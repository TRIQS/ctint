#include "./solver_core.hpp"
#include "./measures.hpp"
#include "./moves/insert.hpp"
#include "./moves/remove.hpp"
#include "./post_process.hpp"
#include "./qmc_config.hpp"
#include "./vertex_factories.hpp"

namespace triqs_ctint {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

    // Set inverse temperature for all $\tau$ points
    tau_t::beta = p.beta;

    // Allocate essential QMC containers
    G0_iw    = block_gf<imfreq>{{p.beta, Fermion, p.n_iw}, p.gf_struct};
    G_iw     = G0_iw;
    Sigma_iw = G0_iw;

    // Allocate containers for dynamical density-density interaction
    if (p.use_D) D0_iw = make_block2_gf<imfreq, matrix_valued>({p.beta, Boson, p.n_iw_dynamical_interactions}, p.gf_struct);

    // Allocate containers for dynamical spin-spin interaction
    if (p.use_Jperp) {
      auto bl  = *(p.gf_struct.begin());
      Jperp_iw = gf<imfreq, matrix_valued>{{p.beta, Boson, p.n_iw_dynamical_interactions}, make_shape(bl.second.size(), bl.second.size())};
    }
  }

  // -------------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params) {

    last_solve_params = solve_params;

    // http://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20ctint
    if (world.rank() == 0)
      std::cout << "\n"
                   "╔╦╗╦═╗╦╔═╗ ╔═╗  ┌─┐┌┬┐┬ ┌┐┌┌┬┐\n"
                   " ║ ╠╦╝║║═╬╗╚═╗  │   │ │ │││ │ \n"
                   " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ ┴ ┘└┘ ┴ \n";

    // Merge constr_params and solve_params
    params_t params(constr_params, solve_params);

    // Prepare shifted non-interacting Green function for QMC
    prepare_G0_shift_tau(params);

    // Reset the containers
    container_set::operator=(container_set{});

    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<mc_weight_t> mc(params.random_name, params.random_seed, params.verbosity);

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
    if (params.measure_average_k) mc.add_measure(measures::average_k{params, qmc_config, &result_set()}, "perturbation order measure");
    if (params.measure_histogram) mc.add_measure(measures::histogram{params, qmc_config, &result_set()}, "perturbation order histogram measure");
    if (params.measure_density) mc.add_measure(measures::density{params, qmc_config, &result_set()}, "density matrix measure");
    if (params.measure_M_tau) mc.add_measure(measures::M_tau{params, qmc_config, &result_set()}, "M_tau measure");
    if (params.measure_M_iw) mc.add_measure(measures::M_iw{params, qmc_config, &result_set()}, "M_iw measure");
    if (params.measure_M4_iw) mc.add_measure(measures::M4_iw{params, qmc_config, &result_set()}, "M4_iw measure");
    if (params.measure_M3pp_iw) mc.add_measure(measures::M3pp_iw{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_iw measure");
    if (params.measure_M3ph_iw) mc.add_measure(measures::M3ph_iw{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_iw measure");
    if (params.measure_M3pp_tau) mc.add_measure(measures::M3pp_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_tau measure");
    if (params.measure_M3ph_tau) mc.add_measure(measures::M3ph_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_tau measure");
    if (params.measure_chi2pp_tau) mc.add_measure(measures::chi2_tau<Chan_t::PP>{params, qmc_config, &result_set()}, "chi2pp_tau measure");
    if (params.measure_chi2ph_tau) mc.add_measure(measures::chi2_tau<Chan_t::PH>{params, qmc_config, &result_set()}, "chi2ph_tau measure");
    if (params.measure_chiAB_tau) mc.add_measure(measures::chiAB_tau{params, qmc_config, &result_set()}, "chiAB_tau measure");

    // Perform QMC run and collect results
    mc.warmup_and_accumulate(params.n_warmup_cycles, params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time));
    mc.collect_results(world);

    if (world.rank() == 0) {
      if (params.measure_average_sign) std::cout << "Average sign: " << average_sign << "\n";
      if (params.measure_average_k) std::cout << "Average perturbation order: " << average_k << "\n";
    }

    // Post Processing
    if (params.post_process) { post_process(params); }
  }

  // -------------------------------------------------------------------------------

  void solver_core::prepare_G0_shift_tau(params_t const &p) {

    // Assert hermiticity of the given Weiss field FIXME
    if (!is_gf_hermitian(G0_iw)) TRIQS_RUNTIME_ERROR << "Please make sure that G0_iw fullfills the hermiticity relation G_ij[iw] = G_ji[-iw]*";

    // Prepare shifted non-interacting Green Function G0_shift_tau for Monte Carlo
    // with renormalization of the chemical potential due to alpha
    g_iw_t G0_inv = inverse(G0_iw);

    // Assert compatibility between gf_struct an alpha
    if (p.gf_struct.size() != p.alpha.size()) TRIQS_RUNTIME_ERROR << "Error: Alpha and gf_struct_t incompatible: Different number of blocks \n";
    for (auto [bl, alpha_bl] : zip(p.gf_struct, p.alpha))
      if (alpha_bl.shape() != make_shape(bl.second.size(), p.n_s)) TRIQS_RUNTIME_ERROR << "Error: Alpha block-shape incompatible with gf_struct \n";

    // Loop over static density-density interaction terms
    for (auto const &term : p.h_int) {

      auto &m = term.monomial;

      if (m[0].indices[0] != m[3].indices[0] or m[1].indices[0] != m[2].indices[0])
        TRIQS_RUNTIME_ERROR << "Interaction term with incompatible block structure: cdag_1 cdag_2 c_2 c_1 required";

      if (!is_densdens_interact(m)) continue;

      auto [bl1, idx1] = get_int_indices(m[0], p.gf_struct);
      auto [bl2, idx2] = get_int_indices(m[1], p.gf_struct);

      // Shift diagonal Green function components according to Eq. (17) of implementation Notes
      dcomplex shift_1 = 0.0;
      dcomplex shift_2 = 0.0;
      for (int s : range(p.n_s)) {
        shift_1 += p.alpha[bl1](idx1, s);
        shift_2 += p.alpha[bl2](idx2, s);
      }
      shift_1 *= dcomplex(term.coef) / p.n_s;
      shift_2 *= dcomplex(term.coef) / p.n_s;

      auto g_1 = slice_target_to_scalar(G0_inv[bl1], idx1, idx1);
      auto g_2 = slice_target_to_scalar(G0_inv[bl2], idx2, idx2);
      g_1(iw_) << g_1(iw_) - shift_2;
      g_2(iw_) << g_2(iw_) - shift_1;
    }

    if (D0_iw) {

      // External loop over blocks
      for (int sig : range(p.n_blocks())) {

        // Get Matrix Rank for block
        int Rank = G0_iw[sig].target_shape()[0];
        for (int i : range(Rank)) {

          // Calculate and subtract term according to notes
          dcomplex term = 0.0;
          for (int sigp : range(p.n_blocks()))
            for (int j : range(Rank))
              for (int s : range(p.n_s)) { term += ((*D0_iw)(sig, sigp)[0](i, j) + (*D0_iw)(sigp, sig)[0](j, i)) * p.alpha[sigp](j, s); }
          auto g = slice_target_to_scalar(G0_inv[sig], i, i);
          g(iw_) << g(iw_) - term / p.n_s;
        }
      }
    }

    // Invert and Fourier transform to imaginary times
    G0_shift_iw = inverse(G0_inv);

    // Known high-frequency moments of G0_shift_iw are assumed to be {0,1}
    auto km = make_zero_tail(G0_shift_iw, 2);
    for (auto &km_bl : km) matrix_view<dcomplex>{km_bl(1, ellipsis())} = 1.0;
    auto tau_mesh = gf_mesh<imtime>{p.beta, Fermion, p.n_tau};
#ifdef GTAU_IS_COMPLEX
    auto [tail, err] = fit_hermitian_tail(G0_shift_iw, km);
    G0_shift_tau     = make_gf_from_fourier(G0_shift_iw, tau_mesh, tail);
#else
    if (!is_gf_real_in_tau(G0_shift_iw, 1e-8)) {
      std::cerr << "WARNING: Assuming real G(tau), but found violation |G(iw) - G*(-iw)| > 1e-8. Making it real in tau.\n";
      G0_shift_iw = make_real_in_tau(G0_shift_iw);
    }
    auto [tail, err] = fit_hermitian_tail(G0_shift_iw, km);
    G0_shift_tau     = real(make_gf_from_fourier(G0_shift_iw, tau_mesh, tail));
#endif
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &p) {

    if (world.rank() == 0)
      std::cout << "\n"
                   "Post-processing ... \n";

    // Calculate M_iw from M_tau (Cast from matrix_real_valued to matrix_valued)
    // Set known_moments to zero, in order to avoid tau-derivative fitting in M_tau
    if (M_tau) {
      M_iw = make_gf_from_fourier(block_gf<imtime, matrix_valued>{*M_tau}, G0_iw[0].mesh(), make_zero_tail(G0_iw));
      M_iw = make_hermitian(M_iw.value());
      for (auto [M_bl, M_hartree_bl] : zip(M_iw.value(), M_hartree.value())) M_bl(iw_) << M_bl[iw_] + M_hartree_bl;
    }

    // Calculate G_iw and Sigma_iw from M_iw
    if (M_iw) {
      G_iw     = G0_shift_iw + G0_shift_iw * M_iw.value() * G0_shift_iw;
      Sigma_iw = inverse(G0_iw) - inverse(G_iw); // Careful, dont use shifted Gf here
    }

    // Calculate M3_iw from M3_tau
    if (M3pp_tau) {
      {
        auto iw_mesh       = gf_mesh<imfreq>{p.beta, Fermion, p.n_iw_M3};
        auto iW_mesh       = gf_mesh<imfreq>{p.beta, Boson, p.n_iW_M3};
        auto iw_mesh_large = gf_mesh<imfreq>{p.beta, Fermion, p.n_iw_M3 + p.n_iW_M3};
        auto M3pp_ferm_iw  = make_gf_from_fourier<0, 1>(M3pp_tau.value(), iw_mesh, iw_mesh_large);
        auto M3pp_del_iW   = make_gf_from_fourier(M3pp_delta.value(), iW_mesh, make_zero_tail(M3pp_delta.value()));
        M3pp_iw            = make_block2_gf(gf_mesh{iw_mesh, iW_mesh}, p.gf_struct);

        // Shift from fermionic to mixed particle-particle frequency notation
        M3pp_iw.value()(bl1_, bl2_)(iw_, iW_)(i_, j_, k_, l_)
           << M3pp_ferm_iw(bl1_, bl2_)(iw_, iW_ - iw_)(i_, j_, k_, l_) + M3pp_del_iW(bl1_, bl2_)(iW_)(i_, j_, k_, l_);

        //// CAUTION! The both times should be fourier transformed with e^{-iwt}
        //// We correct this with an overall minus sign for both frequencies
        //M3pp_iw.value()(bl1_, bl2_)(iw_, iW_)(i_, j_, k_, l_) << M3pp_ferm_iw(bl1_, bl2_)(-iw_, -(iW_ - iw_))(i_, j_, k_, l_) + M3pp_del_iW(bl1_, bl2_)(-iW_)(i_, j_, k_, l_);
      }

      if (M_iw) {
        M2pp_tau       = M2_from_M3<Chan_t::PP>(M3pp_tau.value(), M3pp_delta.value(), M_iw.value(), G0_shift_iw, M_tau.value(), M_hartree.value(),
                                          G0_shift_tau, p.gf_struct);
        chi2pp_new_tau = chi2_from_M2<Chan_t::PP>(M2pp_tau.value(), M_iw.value(), G0_shift_iw, density.value());
        auto iw_mesh   = gf_mesh<imfreq>{p.beta, Boson, p.n_iw_chi2};
        chi2pp_new_iw  = make_gf_from_fourier(chi2pp_new_tau.value(), iw_mesh, make_zero_tail(chi2pp_new_tau.value()));
      }
    }
    if (M3ph_tau) {
      {
        auto iw_mesh       = gf_mesh<imfreq>{p.beta, Fermion, p.n_iw_M3};
        auto iW_mesh       = gf_mesh<imfreq>{p.beta, Boson, p.n_iW_M3};
        auto iw_mesh_large = gf_mesh<imfreq>{p.beta, Fermion, p.n_iw_M3 + p.n_iW_M3};
        auto M3ph_ferm_iw  = make_gf_from_fourier<0, 1>(M3ph_tau.value(), iw_mesh, iw_mesh_large);
        auto M3ph_del_iW   = make_gf_from_fourier(M3ph_delta.value(), iW_mesh, make_zero_tail(M3ph_delta.value()));
        M3ph_iw            = make_block2_gf(gf_mesh{iw_mesh, iW_mesh}, p.gf_struct);

        // Shift from fermionic to mixed particle-hole frequency notation
        // CAUTION! The first time should be fourier transformed with e^{-iwt}
        // We correct this with an overall minus sign for the first frequency
        M3ph_iw.value()(bl1_, bl2_)(iw_, iW_)(i_, j_, k_, l_)
           << M3ph_ferm_iw(bl1_, bl2_)(-iw_, iW_ + iw_)(i_, j_, k_, l_) + M3ph_del_iW(bl1_, bl2_)(iW_)(i_, j_, k_, l_);
      }

      if (M_iw) {
        M2ph_tau       = M2_from_M3<Chan_t::PH>(M3ph_tau.value(), M3ph_delta.value(), M_iw.value(), G0_shift_iw, M_tau.value(), M_hartree.value(),
                                          G0_shift_tau, p.gf_struct);
        chi2ph_new_tau = chi2_from_M2<Chan_t::PH>(M2ph_tau.value(), M_iw.value(), G0_shift_iw, density.value());
        auto iw_mesh   = gf_mesh<imfreq>{p.beta, Boson, p.n_iw_chi2};
        chi2ph_new_iw  = make_gf_from_fourier(chi2ph_new_tau.value(), iw_mesh, make_zero_tail(chi2ph_new_tau.value()));
      }
    }

    // Calculate G2c_iw, F_iw and G2_iw from M4_iw and M_iw
    if (M4_iw and M_iw) G2c_iw = G2c_from_M4(M4_iw.value(), M_iw.value(), G0_shift_iw);
    if (G2c_iw and M_iw) F_iw = F_from_G2c(G2c_iw.value(), G_iw);
    if (G2c_iw and M_iw) G2_iw = G2_from_G2c(G2c_iw.value(), G_iw);

    // Calculate chi3_iw from M3_iw and M_iw
    if (M3pp_iw and M_iw) chi3pp_iw = chi3_from_M3<Chan_t::PP>(M3pp_iw.value(), M_iw.value(), G0_shift_iw, density.value());
    if (M3ph_iw and M_iw) chi3ph_iw = chi3_from_M3<Chan_t::PH>(M3ph_iw.value(), M_iw.value(), G0_shift_iw, density.value());
    if (M3pp_iw_nfft and M_iw) chi3pp_iw_nfft = chi3_from_M3<Chan_t::PP>(M3pp_iw_nfft.value(), M_iw.value(), G0_shift_iw, density.value());
    if (M3ph_iw_nfft and M_iw) chi3ph_iw_nfft = chi3_from_M3<Chan_t::PH>(M3ph_iw_nfft.value(), M_iw.value(), G0_shift_iw, density.value());

    // Calculate chi2_iw from chi2_tau
    auto iw_mesh = gf_mesh<imfreq>{p.beta, Boson, p.n_iw_chi2};
    if (chi2pp_tau) chi2pp_iw = make_gf_from_fourier(chi2pp_tau.value(), iw_mesh, make_zero_tail(chi2pp_tau.value()));
    if (chi2ph_tau) chi2ph_iw = make_gf_from_fourier(chi2ph_tau.value(), iw_mesh, make_zero_tail(chi2ph_tau.value()));

    // Calculate chiAB_iw from chiAB_tau
    if (chiAB_tau) chiAB_iw = make_gf_from_fourier(chiAB_tau.value(), iw_mesh, make_zero_tail(chiAB_tau.value()));
  }

} // namespace triqs_ctint
