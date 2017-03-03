#include "./solver_core.hpp"
#include "./vertex_factories.hpp"
#include "./qmc_config.hpp"
#include "./moves/insert.hpp"
#include "./moves/remove.hpp"
#include "./measures.hpp"
#include "./post_process.hpp"
#include "./fourier_factories.hpp"

namespace triqs_ctint {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

    // Set inverse temperature for all $\tau$ points
    tau_t::beta = p.beta;

    // Allocate essential QMC containers
    G0_iw        = make_block_gf(gf_mesh<imfreq>{p.beta, Fermion, p.n_iw}, p.gf_struct);
    G0_shift_tau = make_block_gf(gf_mesh<imtime>{p.beta, Fermion, p.n_tau}, p.gf_struct);

    // Allocate containers for dynamical density-density interaction // FIXME Block2_gf ?
    if (p.use_D) {
      auto D_block_names = std::vector<std::string>{};
      for (auto const &str1 : p.block_names())
        for (auto const &str2 : p.block_names()) D_block_names.push_back(str1 + "|" + str2);

      std::vector<gf<imfreq, matrix_valued>> v_iw;
      for (auto const &bl1 : p.gf_struct)
        for (auto const &bl2 : p.gf_struct)
          v_iw.emplace_back(
             gf<imfreq, matrix_valued>{{p.beta, Boson, p.n_iw_dynamical_interactions}, make_shape(bl1.second.size(), bl2.second.size())});

      D0_iw = make_block_gf(D_block_names, v_iw);
    }

    // Allocate containers for dynamical spin-spin interaction
    if (p.use_Jperp) {
      auto bl  = *(p.gf_struct.begin());
      Jperp_iw = gf<imfreq, matrix_valued>{{p.beta, Boson, p.n_iw_dynamical_interactions}, make_shape(bl.second.size(), bl.second.size())};
    }
  }

  // -------------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params) {

    if (world.rank() == 0) std::cout << "Welcome to the CT-INT solver" << std::endl << std::endl;

    // Merge constr_params and solve_params
    params_t params(constr_params, solve_params);

    // Prepare shifted non-interacting Green function for QMC
    prepare_G0_shift_tau(params);

    // Reset the containers
    container_set::operator=(container_set{});

    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<double> mc(params.random_name, params.random_seed, 1.0 /* initial sign */, params.verbosity);

    // Capture random number generator
    auto &rng = mc.get_rng();

    // Create Monte-Carlo configuration
    qmc_config_t qmc_config(params, G0_shift_tau);

    // Build vertex factories
    const std::vector<vertex_factory_t> vertex_factories = make_vertex_factories(params, rng, D0_iw, Jperp_iw);

    mc.add_move(moves::insert{&qmc_config, vertex_factories, rng, false}, "insertion");
    mc.add_move(moves::remove{&qmc_config, vertex_factories, rng, false}, "removal");
    if (params.use_double_insertion) {
      mc.add_move(moves::insert{&qmc_config, vertex_factories, rng, true}, "double insertion");
      mc.add_move(moves::remove{&qmc_config, vertex_factories, rng, true}, "double removal");
    }

    // Register all measurements
    if (params.measure_average_sign) mc.add_measure(measures::average_sign{params, qmc_config, &result_set()}, "sign measure");
    if (params.measure_M_tau) mc.add_measure(measures::M_tau{params, qmc_config, &result_set()}, "M_tau measure");
    if (params.measure_M_iw) mc.add_measure(measures::M_iw{params, qmc_config, &result_set()}, "M_iw measure");
    if (params.measure_F_tau) mc.add_measure(measures::F_tau{params, qmc_config, &result_set(), G0_shift_tau}, "F_tau measure");
    if (params.measure_M4_iw) mc.add_measure(measures::M4_iw{params, qmc_config, &result_set()}, "M4_iw measure");
    if (params.measure_M3pp_iw) mc.add_measure(measures::M3_iw<Chan_t::PP>{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_iw measure");
    if (params.measure_M3ph_iw) mc.add_measure(measures::M3_iw<Chan_t::PH>{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_iw measure");
    if (params.measure_M3xph_iw) mc.add_measure(measures::M3_iw<Chan_t::XPH>{params, qmc_config, &result_set(), G0_shift_tau}, "M3xph_iw measure");

    // Perform QMC run and collect results
    mc.warmup_and_accumulate(params.n_warmup_cycles, params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time));
    mc.collect_results(world);

    // Post Processing
    if (params.post_process) post_process(params, qmc_config, &result_set());

    // Write results to file
    triqs::h5::file h5file("ctqmc_out.h5", 'w');
    h5_write(h5file, "", *this);
  }

  // -------------------------------------------------------------------------------

  void solver_core::prepare_G0_shift_tau(params_t const &p) {

    // Prepare shifted non-interacting Green Function G0_shift_tau for Monte Carlo

    auto U = p.get_U();

    // Renormalization of the chemical potential due to alpha
    // Renormalized static interaction by the dynamical one (T.Ayral thesis Eq.11.23a)
    // mu_{sigma i} ---> mu_{sigma i} - (1/n_s) sum_{j sigma' s} ( U_{ij sigma sigma'} - 2.0 int dtau D_{ij sigma sigma'}(tau-0) ) alpha_{s j sigma'}
    // 1/n_s comes immediately from adding an additional sum over s

    decltype(G0_iw) G0_inv = inverse(G0_iw);

    // External loop over blocks
    for (int sig : range(p.n_blocks())) {

      // Get Matrix Rank for block
      int Rank = G0_iw[sig].target_shape()[0];
      for (int i : range(Rank)) {

        // Calculate and subtract term according to equation above
        dcomplex term = 0.0;
        for (int sigp : range(p.n_blocks()))
          for (int j : range(Rank))
            for (int s : range(p.n_s)) {
              term += U(sig, sigp)(i, j) * p.alpha[sigp](j, s);
              if (D0_iw) term += (*D0_iw)[sig * p.n_blocks() + sigp][0](i, j) * p.alpha[sigp](j, s);
            }
        auto g         = slice_target_to_scalar(G0_inv[sig], i, i);
        double h_shift = (p.hartree_shift.size() > 0) ? p.hartree_shift[sig] : 0.0;
        g(iw_) << g(iw_) + h_shift - term / p.n_s;
      }
    }
    // Invert and Fourier transform to imaginary times
    G0_shift_tau = make_gf_from_inverse_fourier(inverse(G0_inv), p.n_tau);
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &p, qmc_config_t const &qmc_config, container_set *results) {

    // Calculate M_iw from M_tau
    if (M_tau) M_iw = make_gf_from_fourier(*M_tau, p.n_iw_M4);
  }

} // namespace triqs_ctint
