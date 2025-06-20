// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

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
      auto [bl, bl_size] = p.gf_struct[0];
      Jperp_iw           = gf<imfreq, matrix_valued>{{p.beta, Boson, p.n_iw_dynamical_interactions}, make_shape(bl_size, bl_size)};
    }
  }

  // -------------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params) {

    last_solve_params = solve_params;

    // Merge constr_params and solve_params
    params_t params(constr_params, solve_params);

    // Open new report stream
    triqs::utility::report_stream report(&std::cout, params.verbosity);

    // http://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20ctint
    report(3) << "\n"
                 "╔╦╗╦═╗╦╔═╗ ╔═╗  ┌─┐┌┬┐┬ ┌┐┌┌┬┐\n"
                 " ║ ╠╦╝║║═╬╗╚═╗  │   │ │ │││ │ \n"
                 " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ ┴ ┘└┘ ┴ \n";

    // Assert hermiticity of the given Weiss field
    if (!is_gf_hermitian(G0_iw)) TRIQS_RUNTIME_ERROR << "Please make sure that G0_iw fullfills the hermiticity relation G_ij[iw] = G_ji[-iw]*";

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

    // Eliminate duplicates from the move types
    if (params.insertion_types.empty()) { params.insertion_types = params.use_double_insertion ? std::vector<int>{1, 2} : std::vector<int>{1}; }
    std::sort(params.insertion_types.begin(), params.insertion_types.end());
    std::ignore = std::unique(params.insertion_types.begin(), params.insertion_types.end());

    for (auto &&num_insert : params.insertion_types) {
      mc.add_move(moves::insert{&qmc_config, vertex_factories, rng, num_insert, params.max_order}, "insert " + std::to_string(num_insert));
      mc.add_move(moves::remove{&qmc_config, vertex_factories, rng, num_insert, params.max_order}, "remove " + std::to_string(num_insert));
    }

    // Register warmup measurements
    mc.add_measure(measures::average_sign{params, qmc_config, &result_set()}, "sign measure", /* enable_timer */ true, /* report */ true);
    mc.add_measure(measures::average_k{params, qmc_config, &result_set()}, "perturbation order measure", /* enable_timer */ true, /* report */ true);

    // Warmup
    report(3) << "\nWarming up ..." << std::endl;
    mc.run(params.n_warmup_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time), /* do_measure */ true);
    double warmup_time_ = mc.get_accumulation_time();

    // Clear warmup measurements
    mc.clear_measures();
    container_set::operator=(container_set{});

    // Register all measurements
    if (params.measure_average_sign)
      mc.add_measure(measures::average_sign{params, qmc_config, &result_set()}, "sign measure", /* enable_timer */ true, /* report */ true);
    if (params.measure_average_k)
      mc.add_measure(measures::average_k{params, qmc_config, &result_set()}, "perturbation order measure", /* enable_timer */ true,
                     /* report */ true);
    if (params.measure_auto_corr_time) mc.add_measure(measures::auto_corr_time{params, qmc_config, &result_set()}, "Auto-correlation time");
    if (params.measure_sign_only) {
      report(3) << "You selected Sign only mode" << std::endl;
    } else {
      if (params.measure_histogram) mc.add_measure(measures::histogram{params, qmc_config, &result_set()}, "perturbation order histogram measure");
      if (params.measure_density) mc.add_measure(measures::density{params, qmc_config, &result_set()}, "density matrix measure");
      if (params.measure_M_tau) mc.add_measure(measures::M_tau{params, qmc_config, &result_set()}, "M_tau measure");
      if (params.measure_M_iw) mc.add_measure(measures::M_iw{params, qmc_config, &result_set()}, "M_iw measure");
      if (params.measure_M4_iw) mc.add_measure(measures::M4_iw{params, qmc_config, &result_set()}, "M4_iw measure");
      if (params.measure_M4pp_iw) mc.add_measure(measures::M4pp_iw{params, qmc_config, &result_set()}, "M4pp_iw measure");
      if (params.measure_M4ph_iw) mc.add_measure(measures::M4ph_iw{params, qmc_config, &result_set()}, "M4ph_iw measure");
      if (params.measure_M3pp_iw) mc.add_measure(measures::M3pp_iw{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_iw measure");
      if (params.measure_M3ph_iw) mc.add_measure(measures::M3ph_iw{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_iw measure");
      if (params.measure_M3pp_tau) mc.add_measure(measures::M3pp_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_tau measure");
      if (params.measure_M3ph_tau) mc.add_measure(measures::M3ph_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_tau measure");
      if (params.measure_M3xph_tau) mc.add_measure(measures::M3xph_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3xph_tau measure");
      if (params.measure_chi2pp_tau) mc.add_measure(measures::chi2_tau<Chan_t::PP>{params, qmc_config, &result_set()}, "chi2pp_tau measure");
      if (params.measure_chi2ph_tau) mc.add_measure(measures::chi2_tau<Chan_t::PH>{params, qmc_config, &result_set()}, "chi2ph_tau measure");
      if (params.measure_chiAB_tau) mc.add_measure(measures::chiAB_tau{params, qmc_config, &result_set()}, "chiAB_tau measure");
    }

    // Perform QMC run and collect results
    report(3) << "\nAccumulating ..." << std::endl;
    mc.run(params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time), /* do_measure */ true);
    double accumulation_time_ = mc.get_accumulation_time();
    mc.collect_results(world);
    warmup_time       = mpi::all_reduce(warmup_time_, world, MPI_MAX);
    accumulation_time = mpi::all_reduce(accumulation_time_, world, MPI_MAX);

    if (params.measure_average_sign) report(3) << "Average sign: " << average_sign << "\n";
    if (params.measure_average_k) report(3) << "Average perturbation order: " << average_k << "\n";
    if (params.measure_auto_corr_time) report(3) << "Auto-correlation time: " << auto_corr_time << "\n";

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
      if (alpha_bl.shape() != make_shape(bl.second, p.n_s)) TRIQS_RUNTIME_ERROR << "Error: Alpha block-shape incompatible with gf_struct \n";

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
    auto tau_mesh = mesh::imtime{p.beta, Fermion, p.n_tau};
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
        auto iW_mesh       = mesh::imfreq{p.beta, Boson, p.n_iW_M3};
        auto iw_mesh       = mesh::imfreq{p.beta, Fermion, p.n_iw_M3};
        auto iw_mesh_large = mesh::imfreq{p.beta, Fermion, p.n_iw_M3 + p.n_iW_M3};
        auto M3pp_ferm_iw  = make_gf_from_fourier<0, 1>(M3pp_tau.value(), iw_mesh, iw_mesh_large);
        auto M3pp_del_iW   = make_gf_from_fourier(M3pp_delta.value(), iW_mesh, make_zero_tail(M3pp_delta.value()));
        M3pp_iw            = make_block2_gf(prod{iW_mesh, iw_mesh}, p.gf_struct);

        // Shift from fermionic to mixed particle-particle frequency notation
        M3pp_iw.value()(bl1_, bl2_)(iW_, iw_)(i_, j_, k_, l_)
           << M3pp_ferm_iw(bl1_, bl2_)(iw_, iW_ - iw_)(i_, j_, k_, l_) + M3pp_del_iW(bl1_, bl2_)(iW_)(i_, j_, k_, l_);

        //// CAUTION! The both times should be fourier transformed with e^{-iwt}
        //// We correct this with an overall minus sign for both frequencies
        //M3pp_iw.value()(bl1_, bl2_)(iW_, iw_)(i_, j_, k_, l_) << M3pp_ferm_iw(bl1_, bl2_)(-iw_, -(iW_ - iw_))(i_, j_, k_, l_) + M3pp_del_iW(bl1_, bl2_)(-iW_)(i_, j_, k_, l_);
      }

      if (M_iw) {
        chi2pp_conn_tau_from_M3 = chi2_conn_from_M3<Chan_t::PP>(M3pp_tau.value(), M3pp_delta.value(), M_iw.value(), G0_shift_iw, M_tau.value(),
                                                                M_hartree.value(), G0_shift_tau);
        chi2pp_tau_from_M3      = chi2_from_chi2_conn<Chan_t::PP>(chi2pp_conn_tau_from_M3.value(), G_iw, density.value());
        auto iw_mesh            = mesh::imfreq{p.beta, Boson, p.n_iw_chi2};
        chi2pp_iw_from_M3       = make_gf_from_fourier(chi2pp_tau_from_M3.value(), iw_mesh, make_zero_tail(chi2pp_tau_from_M3.value()));
      }
    }
    if (M3ph_tau) {
      {
        auto iW_mesh       = mesh::imfreq{p.beta, Boson, p.n_iW_M3};
        auto iw_mesh       = mesh::imfreq{p.beta, Fermion, p.n_iw_M3};
        auto iw_mesh_large = mesh::imfreq{p.beta, Fermion, p.n_iw_M3 + p.n_iW_M3};
        auto M3ph_ferm_iw  = make_gf_from_fourier<0, 1>(M3ph_tau.value(), iw_mesh, iw_mesh_large);
        auto M3ph_del_iW   = make_gf_from_fourier(M3ph_delta.value(), iW_mesh, make_zero_tail(M3ph_delta.value()));
        M3ph_iw            = make_block2_gf(prod{iW_mesh, iw_mesh}, p.gf_struct);

        // Shift from fermionic to mixed particle-hole frequency notation
        // CAUTION! The first time should be fourier transformed with e^{-iwt}
        // We correct this with an overall minus sign for the first frequency
        M3ph_iw.value()(bl1_, bl2_)(iW_, iw_)(i_, j_, k_, l_)
           << M3ph_ferm_iw(bl1_, bl2_)(-iw_, iW_ + iw_)(i_, j_, k_, l_) + M3ph_del_iW(bl1_, bl2_)(iW_)(i_, j_, k_, l_);
      }

      if (M_iw) {
        chi2ph_conn_tau_from_M3 = chi2_conn_from_M3<Chan_t::PH>(M3ph_tau.value(), M3ph_delta.value(), M_iw.value(), G0_shift_iw, M_tau.value(),
                                                                M_hartree.value(), G0_shift_tau);
        chi2ph_tau_from_M3      = chi2_from_chi2_conn<Chan_t::PH>(chi2ph_conn_tau_from_M3.value(), G_iw, density.value());
        auto iw_mesh            = mesh::imfreq{p.beta, Boson, p.n_iw_chi2};
        chi2ph_iw_from_M3       = make_gf_from_fourier(chi2ph_tau_from_M3.value(), iw_mesh, make_zero_tail(chi2ph_tau_from_M3.value()));
      }
    }
    if (M3xph_tau) {
      {
        auto iw_mesh       = mesh::imfreq{p.beta, Fermion, p.n_iw_M3};
        auto iW_mesh       = mesh::imfreq{p.beta, Boson, p.n_iW_M3};
        auto iw_mesh_large = mesh::imfreq{p.beta, Fermion, p.n_iw_M3 + p.n_iW_M3};
        auto M3xph_ferm_iw = make_gf_from_fourier<0, 1>(M3xph_tau.value(), iw_mesh, iw_mesh_large);
        auto M3xph_del_iW  = make_gf_from_fourier(M3xph_delta.value(), iW_mesh, make_zero_tail(M3xph_delta.value()));
        M3xph_iw           = make_block2_gf(mesh::prod{iW_mesh, iw_mesh}, p.gf_struct);

        // Shift from fermionic to mixed particle-hole-cross frequency notation
        // CAUTION! The first time should be fourier transformed with e^{-iwt}
        // We correct this with an overall minus sign for the first frequency
        M3xph_iw.value()(bl1_, bl2_)(iW_, iw_)(i_, j_, k_, l_)
           << M3xph_ferm_iw(bl1_, bl2_)(iW_ + iw_, -iw_)(i_, j_, k_, l_) + M3xph_del_iW(bl1_, bl2_)(iW_)(i_, j_, k_, l_);
      }

      if (M_iw) {
        chi2xph_conn_tau_from_M3 = chi2_conn_from_M3<Chan_t::XPH>(M3xph_tau.value(), M3xph_delta.value(), M_iw.value(), G0_shift_iw, M_tau.value(),
                                                                  M_hartree.value(), G0_shift_tau);
        chi2xph_tau_from_M3      = chi2_from_chi2_conn<Chan_t::XPH>(chi2xph_conn_tau_from_M3.value(), G_iw, density.value());
        auto iw_mesh             = mesh::imfreq{p.beta, Boson, p.n_iw_chi2};
        chi2xph_iw_from_M3       = make_gf_from_fourier(chi2xph_tau_from_M3.value(), iw_mesh, make_zero_tail(chi2xph_tau_from_M3.value()));
      }
    }

    // Calculate G2_conn_iw, F_iw and G2_iw from M4_iw and M_iw
    if (M4_iw and M_iw) G2_conn_iw = G2_conn_from_M4(M4_iw.value(), M_iw.value(), G0_shift_iw);
    if (M4pp_iw and M_iw) G2pp_conn_iw = G2pp_conn_from_M4pp(M4pp_iw.value(), M_iw.value(), G0_shift_iw);
    if (M4ph_iw and M_iw) G2ph_conn_iw = G2ph_conn_from_M4ph(M4ph_iw.value(), M_iw.value(), G0_shift_iw);

    if (G2_conn_iw and M_iw) F_iw = F_from_G2c(G2_conn_iw.value(), G_iw);
    if (G2pp_conn_iw and M_iw) Fpp_iw = Fpp_from_G2pp_conn(G2pp_conn_iw.value(), G_iw);
    if (G2ph_conn_iw and M_iw) Fph_iw = Fph_from_G2ph_conn(G2ph_conn_iw.value(), G_iw);

    if (G2_conn_iw and M_iw) G2_iw = G2_from_G2c(G2_conn_iw.value(), G_iw);
    if (G2pp_conn_iw and M_iw) G2pp_iw = G2pp_from_G2pp_conn(G2pp_conn_iw.value(), G_iw);
    if (G2ph_conn_iw and M_iw) G2ph_iw = G2ph_from_G2ph_conn(G2ph_conn_iw.value(), G_iw);

    // Calculate chi3_iw from M3_iw and M_iw
    if (M3pp_iw and M_iw) chi3pp_iw = chi3_from_M3<Chan_t::PP>(M3pp_iw.value(), M_iw.value(), G0_shift_iw, density.value(), M_hartree.value());
    if (M3ph_iw and M_iw) chi3ph_iw = chi3_from_M3<Chan_t::PH>(M3ph_iw.value(), M_iw.value(), G0_shift_iw, density.value(), M_hartree.value());
    if (M3xph_iw and M_iw) chi3xph_iw = chi3_from_M3<Chan_t::XPH>(M3xph_iw.value(), M_iw.value(), G0_shift_iw, density.value(), M_hartree.value());
    if (M3pp_iw_nfft and M_iw)
      chi3pp_iw_nfft = chi3_from_M3<Chan_t::PP>(M3pp_iw_nfft.value(), M_iw.value(), G0_shift_iw, density.value(), M_hartree.value());
    if (M3ph_iw_nfft and M_iw)
      chi3ph_iw_nfft = chi3_from_M3<Chan_t::PH>(M3ph_iw_nfft.value(), M_iw.value(), G0_shift_iw, density.value(), M_hartree.value());

    // Calculate chi2_iw from chi2_tau
    auto iw_mesh = mesh::imfreq{p.beta, Boson, p.n_iw_chi2};
    if (chi2pp_tau) chi2pp_iw = make_gf_from_fourier(chi2pp_tau.value(), iw_mesh, make_zero_tail(chi2pp_tau.value()));
    if (chi2ph_tau) chi2ph_iw = make_gf_from_fourier(chi2ph_tau.value(), iw_mesh, make_zero_tail(chi2ph_tau.value()));

    // Calculate chiAB_iw from chiAB_tau
    if (chiAB_tau) chiAB_iw = make_gf_from_fourier(chiAB_tau.value(), iw_mesh, make_zero_tail(chiAB_tau.value()));
  }

} // namespace triqs_ctint
