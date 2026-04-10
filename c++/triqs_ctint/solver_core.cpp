// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./solver_core.hpp"
#include "./measures.hpp"

#include <chrono>
#include "./moves/insert.hpp"
#include "./moves/remove.hpp"
#include "./moves/spinflip.hpp"
#include "./post_process.hpp"
#include "./qmc_config.hpp"
#include "./vertex_factories.hpp"

namespace triqs_ctint {

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

    // Set inverse temperature for all $\tau$ points
    tau_t::set_beta(p.beta);

    // Allocate essential QMC containers on DLR mesh (symmetrize=true for hermiticity checks)
    G0_iw        = g_iw_t{{p.beta, Fermion, p.dlr_wmax, p.dlr_eps, true}, p.gf_struct};
    G0_iw_inv    = G0_iw;
    G_iw         = G0_iw;
    Sigma_dyn_iw = G0_iw;

    // Allocate containers for dynamical density-density interaction (DLR bosonic mesh)
    if (p.use_D) {
      auto D_mesh = mesh::dlr_imfreq{p.beta, Boson, p.dlr_wmax, p.dlr_eps, true};
      auto bl_sz  = p.gf_struct[0].second;
      auto bl     = p.block_names();
      auto g0     = gf<mesh::dlr_imfreq, matrix_valued>{D_mesh, make_shape(bl_sz, bl_sz)};
      D0_iw       = block2_gf<mesh::dlr_imfreq, matrix_valued>{{bl, bl}, std::vector(bl.size(), std::vector(bl.size(), g0))};
    }

    // Allocate containers for dynamical spin-spin interaction (DLR bosonic mesh)
    if (p.use_Jperp) {
      auto [bl, bl_size] = p.gf_struct[0];
      Jperp_iw = gf<mesh::dlr_imfreq, matrix_valued>{{p.beta, Boson, p.dlr_wmax, p.dlr_eps, true}, make_shape(bl_size, bl_size)};
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

    // Calculate the Weiss field inverse
    G0_iw_inv = inverse(G0_iw);

    // Prepare shifted non-interacting Green function for QMC
    prepare_G0_shift_iw(params);

    // --- Calculate G0_shift_tau from G0_shift_iw via DLR
    auto G0_shift_dlr = make_gf_dlr(G0_shift_iw);
#ifdef GTAU_IS_COMPLEX
    G0_shift_tau = make_gf_imtime(G0_shift_dlr, params.n_tau);
#else
    G0_shift_tau = real(make_gf_imtime(G0_shift_dlr, params.n_tau));
#endif

    // Reset the containers
    container_set::operator=(container_set{});

    // Construct the generic Monte-Carlo solver
    triqs::mc_tools::mc_generic<mc_weight_t> mc(params.random_name, params.random_seed, world, params.verbosity);

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
    if (params.use_auxiliary_spin_flip) {
      TRIQS_ASSERT2(params.n_s == 2, "ERROR: Auxiliary spin-flip move requires n_s = 2.");
      mc.add_move(moves::spinflip{&qmc_config, vertex_factories, rng, /* n_spinflips */ 1, params.max_order}, "auxiliary spin-flip");
    }

    // Determine cycle length for warmup phase
    int warmup_cycle_length = (params.length_cycle >= 0) ? params.length_cycle : 100;
    bool auto_warmup        = (params.n_warmup_cycles < 0);

    // Warmup
    if (auto_warmup) {
      report(3) << "\nWarming up (auto) ..." << std::endl;
      mc.set_verbosity(0);

      // Automatic warmup: run until perturbation order stabilizes
      double mean_k       = 0.0;
      double prev_mean_k  = 0.0;
      mc_weight_t sign_sum = 0.0;
      int64_t n_acc       = 0;
      int n_stable        = 0;
      bool converged      = false;
      double next_print   = 2.0; // first periodic print after 2 seconds

      constexpr int check_interval    = 100;
      constexpr int min_warmup        = 100;
      constexpr double rtol           = 0.03;
      constexpr int n_stable_required = 3;

      auto clock_cb = triqs::utility::clock_callback(params.max_time);
      auto t0       = std::chrono::steady_clock::now();

      auto after_duty = [&]() {
        double k = static_cast<double>(qmc_config.perturbation_order());
        ++n_acc;
        mean_k += (k - mean_k) / n_acc; // online mean
        sign_sum += mc.get_sign();
      };

      auto stop_cb = [&]() -> bool {
        if (clock_cb()) return true;
        if (n_acc < min_warmup || n_acc % check_interval != 0) return false;

        // Periodic status print
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (elapsed > next_print) {
          next_print = 1.25 * elapsed + 2.0;
          report(3) << "  mean k = " << mean_k << ", sign = " << std::real(sign_sum / static_cast<double>(n_acc)) << " (" << n_acc << " cycles)\n";
        }

        // Compute relative change of running mean
        double rel_change = std::abs(mean_k - prev_mean_k) / std::max(std::abs(mean_k), 1.0);
        prev_mean_k       = mean_k;

        if (rel_change < rtol)
          ++n_stable;
        else
          n_stable = 0;

        bool local_converged  = (n_stable >= n_stable_required);
        int global_converged  = mpi::all_reduce(static_cast<int>(local_converged), world, MPI_MIN);
        converged             = (global_converged == 1);
        return converged;
      };

      typename decltype(mc)::run_param_t rp;
      rp.ncycles             = params.max_warmup_cycles;
      rp.cycle_length        = warmup_cycle_length;
      rp.stop_callback       = stop_cb;
      rp.after_cycle_duty    = after_duty;
      rp.comm                = world;
      rp.enable_measures     = false;
      rp.propagate_exception = params.rethrow_exception;
      mc.run(rp);

      if (!converged)
        report(1) << "  WARNING: did not converge after " << params.max_warmup_cycles << " cycles\n";
      else
        report(3) << "  mean k = " << mean_k << " -> converged\n";

      mc.set_verbosity(params.verbosity);
    } else {
      report(3) << "\nWarming up ..." << std::endl;
      // Register warmup measurements for fixed warmup progress display
      mc.add_measure(measures::average_sign{params, qmc_config, &result_set()}, "sign measure", /* enable_timer */ true, /* report */ true);
      mc.add_measure(measures::average_k{params, qmc_config, &result_set()}, "perturbation order measure", /* enable_timer */ true, /* report */ true);
      // Fixed warmup: user-specified number of cycles
      mc.run(params.n_warmup_cycles, warmup_cycle_length, triqs::utility::clock_callback(params.max_time), /* do_measure */ true);
      mc.clear_measures();
    }
    auto warmup_cycles_done_ = mc.get_current_cycle_number();
    double warmup_time_      = mc.get_accumulation_time();

    container_set::operator=(container_set{});
    warmup_cycles_done = warmup_cycles_done_;

    // Automatic length_cycle calibration
    int effective_length_cycle = params.length_cycle;
    bool auto_length_cycle     = (params.length_cycle < 0);
    if (auto_length_cycle) {
      report(3) << "\nCalibrating length_cycle ..." << std::endl;
      mc.set_verbosity(0);

      // Register densities measure for calibration (gives proper auto_corr_time including density correlations)
      mc.add_measure(measures::densities{params, qmc_config, &result_set()}, "calibration densities");

      // Log-binning on perturbation order for convergence detection only
      triqs::stat::log_binning<dcomplex> k_acc(dcomplex{0.0}, -1);
      int64_t calib_count = 0;
      double prev_tau     = -1.0;
      int n_stable        = 0;
      double next_print   = 2.0; // first periodic print after 2 seconds

      constexpr int check_interval    = 500;
      constexpr int min_calib         = 1000;
      constexpr int max_calib         = 100000;
      constexpr double tau_rtol       = 0.1;
      constexpr int n_stable_required = 3;

      auto clock_cb = triqs::utility::clock_callback(params.max_time);
      auto t0       = std::chrono::steady_clock::now();

      auto after_duty = [&]() {
        k_acc << dcomplex(qmc_config.perturbation_order());
        ++calib_count;
      };

      auto stop_cb = [&]() -> bool {
        if (clock_cb()) return true;
        if (calib_count < min_calib || calib_count % check_interval != 0) return false;

        auto [mean, errs, taus, effs] = k_acc.mean_errors_and_taus(world);
        double tau = taus.empty() ? 0.0 : std::real(taus.back());

        // Periodic status print
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (elapsed > next_print) {
          next_print = 1.25 * elapsed + 2.0;
          report(3) << "  tau_ac = " << tau << " (" << calib_count << " cycles)\n";
        }

        if (prev_tau >= 0) {
          double rel_change = std::abs(tau - prev_tau) / std::max(tau, 1.0);
          if (rel_change < tau_rtol)
            ++n_stable;
          else
            n_stable = 0;
        }
        prev_tau = tau;
        return (n_stable >= n_stable_required);
      };

      typename decltype(mc)::run_param_t rp;
      rp.ncycles             = max_calib;
      rp.cycle_length        = 1;
      rp.stop_callback       = stop_cb;
      rp.after_cycle_duty    = after_duty;
      rp.comm                = world;
      rp.enable_measures     = true;
      rp.propagate_exception = params.rethrow_exception;
      mc.run(rp);

      // Collect results to compute auto_corr_time from the densities measure
      mc.collect_results(world);
      double tau_raw = auto_corr_time; // from container_set, includes density + perturbation order correlations

      // Set length_cycle so that effective autocorrelation ≈ target_auto_corr_time
      effective_length_cycle = std::max(1, static_cast<int>(std::ceil(tau_raw / params.target_auto_corr_time)));
      effective_length_cycle = std::min(effective_length_cycle, params.max_length_cycle);

      report(3) << "  tau_ac = " << tau_raw << " -> length_cycle = " << effective_length_cycle << "\n";

      // Clean up calibration phase
      mc.clear_measures();
      auto saved = warmup_cycles_done;
      container_set::operator=(container_set{});
      warmup_cycles_done = saved;
      mc.set_verbosity(params.verbosity);
      report(3) << std::endl;
    }
    length_cycle_used = effective_length_cycle;

    // Auto-enable density_matrix measurement when any M3 measurement is active (needed for chi3 post-processing)
    if (params.measure_M3pp_iw || params.measure_M3ph_iw || params.measure_M3pp_iw_full || params.measure_M3ph_iw_full
        || params.measure_M3pp_tau || params.measure_M3ph_tau || params.measure_M3xph_tau)
      params.measure_density_matrix = true;

    // Register all measurements
    if (params.measure_average_sign)
      mc.add_measure(measures::average_sign{params, qmc_config, &result_set()}, "sign measure", /* enable_timer */ true, /* report */ true);
    if (params.measure_average_k)
      mc.add_measure(measures::average_k{params, qmc_config, &result_set()}, "perturbation order measure", /* enable_timer */ true,
                     /* report */ true);
    mc.add_measure(measures::densities{params, qmc_config, &result_set()}, "densities measure"); // always active for auto-correlation time
    if (params.measure_sign_only) {
      report(3) << "You selected Sign only mode" << std::endl;
    } else {
      if (params.measure_histogram) mc.add_measure(measures::histogram{params, qmc_config, &result_set()}, "perturbation order histogram measure");
      if (params.measure_density_matrix) mc.add_measure(measures::density_matrix{params, qmc_config, &result_set()}, "density matrix measure");
      if (params.measure_M_tau) mc.add_measure(measures::M_tau{params, qmc_config, &result_set()}, "M_tau measure");
      if (params.measure_M_iw) mc.add_measure(measures::M_iw{params, qmc_config, &result_set()}, "M_iw measure");
      if (params.measure_M4_iw) mc.add_measure(measures::M4_iw{params, qmc_config, &result_set()}, "M4_iw measure");
      if (params.measure_M4pp_iw) mc.add_measure(measures::M4pp_iw{params, qmc_config, &result_set()}, "M4pp_iw measure");
      if (params.measure_M4ph_iw) mc.add_measure(measures::M4ph_iw{params, qmc_config, &result_set()}, "M4ph_iw measure");
      if (params.measure_M3pp_iw || params.measure_M3ph_iw || params.measure_M3pp_iw_full || params.measure_M3ph_iw_full)
        report(3) << "Constructing DLR2D mesh (beta = " << params.beta << ", dlr_wmax = " << params.dlr_wmax << ", dlr_eps = " << params.dlr_eps
                  << ") ..." << std::endl;
      if (params.measure_M3pp_iw) mc.add_measure(measures::M3pp_iw{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_iw measure");
      if (params.measure_M3ph_iw) mc.add_measure(measures::M3ph_iw{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_iw measure");
      if (params.measure_M3pp_iw_full) mc.add_measure(measures::M3pp_iw_full{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_iw_full measure");
      if (params.measure_M3ph_iw_full) mc.add_measure(measures::M3ph_iw_full{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_iw_full measure");
      if (params.measure_M3pp_tau) mc.add_measure(measures::M3pp_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3pp_tau measure");
      if (params.measure_M3ph_tau) mc.add_measure(measures::M3ph_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3ph_tau measure");
      if (params.measure_M3xph_tau) mc.add_measure(measures::M3xph_tau{params, qmc_config, &result_set(), G0_shift_tau}, "M3xph_tau measure");
      if (params.measure_chi2pp_tau) mc.add_measure(measures::chi2_tau<Chan_t::PP>{params, qmc_config, &result_set()}, "chi2pp_tau measure");
      if (params.measure_chi2ph_tau) mc.add_measure(measures::chi2_tau<Chan_t::PH>{params, qmc_config, &result_set()}, "chi2ph_tau measure");
      if (params.measure_chiAB_tau) mc.add_measure(measures::chiAB_tau{params, qmc_config, &result_set()}, "chiAB_tau measure");
      if (params.measure_static_obs) mc.add_measure(measures::static_obs{params, qmc_config, &result_set()}, "static_obs measure");
    }

    // Perform QMC run and collect results
    report(3) << "\nAccumulating ..." << std::endl;
    mc.run(params.n_cycles, effective_length_cycle, triqs::utility::clock_callback(params.max_time), /* do_measure */ true);
    double accumulation_time_ = mc.get_accumulation_time();
    mc.collect_results(world);
    warmup_time       = mpi::all_reduce(warmup_time_, world, MPI_MAX);
    accumulation_time = mpi::all_reduce(accumulation_time_, world, MPI_MAX);

    if (params.measure_average_sign) report(3) << "Average sign: " << average_sign << "\n";
    if (params.measure_average_k) report(3) << "Average perturbation order: " << average_k << "\n";
    report(3) << "Auto-correlation time: " << auto_corr_time << "\n";
    report(3) << "Warmup cycles: " << warmup_cycles_done << (auto_warmup ? " (auto)" : "") << "\n";
    report(3) << "Length cycle: " << length_cycle_used << (auto_length_cycle ? " (auto)" : "") << "\n";

    // Post Processing
    if (params.post_process) { post_process(params); }
  }

  // -------------------------------------------------------------------------------

  // Prepare shifted non-interacting Green Function G0_shift_iw for Monte Carlo
  // with renormalization of the chemical potential due to alpha
  void solver_core::prepare_G0_shift_iw(params_t const &p) {

    // Container that will hold the inverse of the shifted Green function
    g_iw_t G0_shift_iw_inv = G0_iw_inv;

    // Assert compatibility between alpha tensor and h_int (+ D0 extension)
    long n_h_int = std::distance(p.h_int.begin(), p.h_int.end());
    long n_D0_total = 0;
    if (D0_iw) {
      long n_bl = p.n_blocks();
      long R    = (*D0_iw)(0, 0).target_shape()[0];
      n_D0_total = n_bl * n_bl * R * R;
    }
    if (p.alpha.shape() != std::array<long, 4>{n_h_int + n_D0_total, 2, 2, p.n_s})
      TRIQS_RUNTIME_ERROR << "Alpha tensor shape " << p.alpha.shape() << " incompatible with h_int (" << n_h_int << " terms) + D0 (" << n_D0_total
                          << " entries)\n";

    // Loop over static density-density interaction terms
    for (auto const &[n, term] : enumerate(p.h_int)) {
      auto &m = term.monomial;

      if (m[0].indices[0] != m[3].indices[0] or m[1].indices[0] != m[2].indices[0])
        TRIQS_RUNTIME_ERROR << "Interaction term with incompatible block structure: cdag_1 cdag_2 c_2 c_1 required";

      auto [bl_0, idx_cdag_0] = get_int_indices(m[0], p.gf_struct);
      auto [bl_1, idx_cdag_1] = get_int_indices(m[1], p.gf_struct);
      auto [bl_c_1, idx_c_1]  = get_int_indices(m[2], p.gf_struct);
      auto [bl_c_0, idx_c_0]  = get_int_indices(m[3], p.gf_struct);

      // Shift equal-time Green function components according to alpha-tensor
      for (long s : range(p.n_s)) {
        G0_shift_iw_inv[bl_0].data()(range::all, idx_cdag_0, idx_c_0) -= p.alpha(n, 1, 1, s) * U_scalar_t(term.coef) / p.n_s;
        G0_shift_iw_inv[bl_1].data()(range::all, idx_cdag_1, idx_c_1) -= p.alpha(n, 0, 0, s) * U_scalar_t(term.coef) / p.n_s;
        if (bl_0 == bl_1) { // Cross terms are only possible for equal blocks
          G0_shift_iw_inv[bl_0].data()(range::all, idx_cdag_0, idx_c_1) += p.alpha(n, 1, 0, s) * U_scalar_t(term.coef) / p.n_s;
          G0_shift_iw_inv[bl_1].data()(range::all, idx_cdag_1, idx_c_0) += p.alpha(n, 0, 1, s) * U_scalar_t(term.coef) / p.n_s;
        }
      }
    }

    if (D0_iw) {

      int n_bl = p.n_blocks();
      int R    = (*D0_iw)(0, 0).target_shape()[0];

      // Extract per-orbital density from D0 alpha entries (averaged over aux spin s)
      // D0 alpha label = n_h_int + sigp * n_bl * R * R + bl2 * R * R + j * R + b
      // alpha[label, 0, 0, s] = density[sigp][j,j] + delta_shift(s)
      // Averaging over s cancels the delta shift
      std::vector<std::vector<g_tau_scalar_t>> block_density(n_bl);
      for (int sigp : range(n_bl)) {
        block_density[sigp].resize(R, g_tau_scalar_t{0});
        for (int j : range(R)) {
          long label = n_h_int + sigp * n_bl * R * R + 0 * R * R + j * R + 0;
          for (int s : range(p.n_s)) block_density[sigp][j] += p.alpha(label, 0, 0, s);
          block_density[sigp][j] /= p.n_s;
        }
      }

      // Precompute D0 at iw=0 (static Hartree contribution) for all block pairs
      std::vector<std::vector<matrix<dcomplex>>> D0_static(n_bl, std::vector<matrix<dcomplex>>(n_bl));
      for (int b1 : range(n_bl))
        for (int b2 : range(n_bl)) D0_static[b1][b2] = make_gf_imfreq(make_gf_dlr((*D0_iw)(b1, b2)), 1)(0);

      for (int sig : range(n_bl)) {
        for (int i : range(R)) {
          dcomplex shift = 0.0;
          for (int sigp : range(n_bl))
            for (int j : range(R)) shift += (D0_static[sig][sigp](i, j) + D0_static[sigp][sig](j, i)) * block_density[sigp][j];
          G0_shift_iw_inv[sig].data()(range::all, i, i) -= shift;
        }
      }
    }

    // Invert
    G0_shift_iw = inverse(G0_shift_iw_inv);
  }

  // -------------------------------------------------------------------------------

  void solver_core::post_process(params_t const &p) {

    if (world.rank() == 0)
      std::cout << "\n"
                   "Post-processing ... \n";

    // --- Determine M_iw (M_dyn only, no Hartree) from either NFFT or FT of M_tau
    if (M_iw_nfft) {
      M_iw = g_iw_t{M_iw_nfft.value()}; // Direct DLR measurement (M_dyn only)
    } else if (M_tau) {
      // DLR fit of M_tau gives M_dyn only (continuous part, no equal-time delta)
      // Cast from g_tau_t::target_t (possibly matrix_real_valued) to matrix_valued for DLR fit
      auto const &mt = *M_tau;
      std::vector<gf<imtime, matrix_valued>> gf_vec;
      for (int b = 0; b < mt.size(); ++b) gf_vec.emplace_back(mt[b]);
      auto M_tau_cast = block_gf<imtime, matrix_valued>{mt.block_names(), std::move(gf_vec)};
      auto M_dlr      = fit_gf_dlr(M_tau_cast, p.dlr_wmax, p.dlr_eps, true);
      M_iw            = g_iw_t{make_gf_dlr_imfreq(M_dlr)};
    }

    // Helper: build M_full = M_dyn + M_hartree (add constant Hartree term to each DLR point)
    auto make_M_full = [&]() {
      g_iw_t M_full = M_iw.value();
      for (auto [M_bl, M_h_bl] : zip(M_full, M_hartree.value()))
        for (long d = 0; d < M_bl.mesh().size(); ++d) M_bl.data()(d, ellipsis()) += M_h_bl;
      return M_full;
    };

    // --- Calculate G_iw and Sigma decomposition from M_iw
    if (M_iw) {
      // Re-initialize G_iw and Sigma_dyn_iw (were reset by container_set::operator=)
      G_iw         = G0_shift_iw;
      Sigma_dyn_iw = G0_shift_iw;

      auto M_full = make_M_full();

      // G_iw = G0_shift + G0_shift * M_full * G0_shift (Dyson equation at DLR points)
      for (auto [G_bl, G0_bl, M_bl] : zip(G_iw, G0_shift_iw, M_full))
        for (auto iw : G_bl.mesh()) G_bl[iw] = G0_bl[iw] + G0_bl[iw] * M_bl[iw] * G0_bl[iw];

      // Sigma_dyn = G0_shift^{-1} - G^{-1} - M_hartree (decays to zero)
      auto G0_shift_iw_inv = inverse(G0_shift_iw);
      auto G_iw_inv        = inverse(G_iw);
      for (auto [S_bl, G0s_bl, Gi_bl, M_h_bl] : zip(Sigma_dyn_iw, G0_shift_iw_inv, G_iw_inv, M_hartree.value()))
        for (auto iw : S_bl.mesh()) S_bl[iw] = G0s_bl[iw] - Gi_bl[iw] - M_h_bl;

      // Sigma_hartree = Sigma_alpha + M_hartree = (G0^{-1} - G0_shift^{-1}) + M_hartree
      Sigma_hartree = make_block_vector<g_tau_scalar_t>(p.gf_struct);
      for (int bl : range(p.n_blocks())) {
        auto alpha_shift = matrix<dcomplex>(G0_iw_inv[bl].data()(0, ellipsis()) - G0_shift_iw_inv[bl].data()(0, ellipsis()));
#ifdef GTAU_IS_COMPLEX
        Sigma_hartree.value()[bl] = matrix<g_tau_scalar_t>(alpha_shift) + M_hartree.value()[bl];
#else
        Sigma_hartree.value()[bl] = real(alpha_shift) + M_hartree.value()[bl];
#endif
      }
    }

    // --- Convert DLR quantities to regular imfreq for higher-order post-processing
    int n_iw_pp = std::max({p.n_iw_M4 + p.n_iW_M4, p.n_iw_M3 + p.n_iW_M3}) + 10;
    // DLR2D NFFT measurements may require a larger mesh to cover all DLR2D frequencies
    if (M3pp_iw_nfft) n_iw_pp = std::max(n_iw_pp, static_cast<int>(M3pp_iw_nfft.value()(0, 0).mesh().max_n()) + 1);
    if (M3ph_iw_nfft) n_iw_pp = std::max(n_iw_pp, static_cast<int>(M3ph_iw_nfft.value()(0, 0).mesh().max_n()) + 1);
    g_reg_iw_t G0_shift_iw_reg = make_gf_imfreq(G0_shift_iw, n_iw_pp);

    std::optional<g_reg_iw_t> M_iw_reg;
    if (M_iw) M_iw_reg = make_gf_imfreq(make_M_full(), n_iw_pp);
    g_reg_iw_t G_iw_reg = make_gf_imfreq(G_iw, n_iw_pp);

    // Calculate M3_iw from M3_tau
    if (M3pp_tau) {
      {
        auto iw_mesh       = mesh::imfreq{p.beta, Fermion, p.n_iw_M3};
        auto M3pp_ferm_iw  = make_gf_from_fourier<0, 1>(M3pp_tau.value(), iw_mesh, iw_mesh);
        auto M3pp_del_iW   = make_gf_from_fourier(M3pp_delta.value(), mesh::imfreq{p.beta, Boson, p.n_iW_M3}, make_zero_tail(M3pp_delta.value()));
        M3pp_iw            = make_block2_gf(prod{iw_mesh, iw_mesh}, p.gf_struct);

        // Convert to fermionic frequency notation (iw1, iw2)
        // CAUTION! Both times should be Fourier transformed with e^{-iwt}
        // We correct this with an overall minus sign for both frequencies
        // The delta term depends on the bosonic transfer frequency Omega = iw1 + iw2
        M3pp_iw.value()(bl1_, bl2_)(iw1_, iw2_)(i_, j_, k_, l_)
           << M3pp_ferm_iw(bl1_, bl2_)(-iw1_, -iw2_)(i_, j_, k_, l_) + M3pp_del_iW(bl1_, bl2_)(-(iw1_ + iw2_))(i_, j_, k_, l_);
      }

    }
    if (M3ph_tau) {
      {
        auto iw_mesh       = mesh::imfreq{p.beta, Fermion, p.n_iw_M3};
        auto M3ph_ferm_iw  = make_gf_from_fourier<0, 1>(M3ph_tau.value(), iw_mesh, iw_mesh);
        auto M3ph_del_iW   = make_gf_from_fourier(M3ph_delta.value(), mesh::imfreq{p.beta, Boson, p.n_iW_M3}, make_zero_tail(M3ph_delta.value()));
        M3ph_iw            = make_block2_gf(prod{iw_mesh, iw_mesh}, p.gf_struct);

        // Convert to fermionic frequency notation (iw1, iw2)
        // CAUTION! The first time should be Fourier transformed with e^{-iwt}
        // We correct this with an overall minus sign for the first frequency
        // The delta term depends on the bosonic transfer frequency Omega = iw2 - iw1
        M3ph_iw.value()(bl1_, bl2_)(iw1_, iw2_)(i_, j_, k_, l_)
           << M3ph_ferm_iw(bl1_, bl2_)(-iw1_, iw2_)(i_, j_, k_, l_) + M3ph_del_iW(bl1_, bl2_)(iw2_ - iw1_)(i_, j_, k_, l_);
      }

    }
    if (M3xph_tau) {
      {
        auto iw_mesh        = mesh::imfreq{p.beta, Fermion, p.n_iw_M3};
        auto M3xph_ferm_iw  = make_gf_from_fourier<0, 1>(M3xph_tau.value(), iw_mesh, iw_mesh);
        auto M3xph_del_iW   = make_gf_from_fourier(M3xph_delta.value(), mesh::imfreq{p.beta, Boson, p.n_iW_M3}, make_zero_tail(M3xph_delta.value()));
        M3xph_iw            = make_block2_gf(mesh::prod{iw_mesh, iw_mesh}, p.gf_struct);

        // Convert to fermionic frequency notation (iw1, iw2)
        // CAUTION! The first time should be Fourier transformed with e^{-iwt}
        // We correct this with an overall minus sign for the first frequency
        // The delta term depends on the bosonic transfer frequency Omega = iw2 - iw1
        M3xph_iw.value()(bl1_, bl2_)(iw1_, iw2_)(i_, j_, k_, l_)
           << M3xph_ferm_iw(bl1_, bl2_)(iw2_, -iw1_)(i_, j_, k_, l_) + M3xph_del_iW(bl1_, bl2_)(iw2_ - iw1_)(i_, j_, k_, l_);
      }

    }

    // Calculate G2_conn_iw, F_iw and G2_iw from M4_iw and M_iw (using regular imfreq quantities)
    if (M4_iw and M_iw_reg) G2_conn_iw = G2_conn_from_M4(M4_iw.value(), M_iw_reg.value(), G0_shift_iw_reg);
    if (M4pp_iw and M_iw_reg) G2pp_conn_iw = G2pp_conn_from_M4pp(M4pp_iw.value(), M_iw_reg.value(), G0_shift_iw_reg);
    if (M4ph_iw and M_iw_reg) G2ph_conn_iw = G2ph_conn_from_M4ph(M4ph_iw.value(), M_iw_reg.value(), G0_shift_iw_reg);

    if (G2_conn_iw and M_iw_reg) F_iw = F_from_G2c(G2_conn_iw.value(), G_iw_reg);
    if (G2pp_conn_iw and M_iw_reg) Fpp_iw = Fpp_from_G2pp_conn(G2pp_conn_iw.value(), G_iw_reg);
    if (G2ph_conn_iw and M_iw_reg) Fph_iw = Fph_from_G2ph_conn(G2ph_conn_iw.value(), G_iw_reg);

    if (G2_conn_iw and M_iw_reg) G2_iw = G2_from_G2c(G2_conn_iw.value(), G_iw_reg);
    if (G2pp_conn_iw and M_iw_reg) G2pp_iw = G2pp_from_G2pp_conn(G2pp_conn_iw.value(), G_iw_reg);
    if (G2ph_conn_iw and M_iw_reg) G2ph_iw = G2ph_from_G2ph_conn(G2ph_conn_iw.value(), G_iw_reg);

    // Calculate chi3_iw from M3_iw and M_iw (using regular imfreq quantities)
    if (M3pp_iw and M_iw_reg and density_matrix) chi3pp_iw = chi3_from_M3<Chan_t::PP>(M3pp_iw.value(), M_iw_reg.value(), G0_shift_iw_reg, density_matrix.value(), M_hartree.value());
    if (M3ph_iw and M_iw_reg and density_matrix) chi3ph_iw = chi3_from_M3<Chan_t::PH>(M3ph_iw.value(), M_iw_reg.value(), G0_shift_iw_reg, density_matrix.value(), M_hartree.value());
    if (M3xph_iw and M_iw_reg and density_matrix) chi3xph_iw = chi3_from_M3<Chan_t::XPH>(M3xph_iw.value(), M_iw_reg.value(), G0_shift_iw_reg, density_matrix.value(), M_hartree.value());
    if (M3pp_iw_nfft and M_iw_reg and density_matrix)
      chi3pp_iw_nfft = chi3_from_M3<Chan_t::PP>(M3pp_iw_nfft.value(), M_iw_reg.value(), G0_shift_iw_reg, density_matrix.value(), M_hartree.value());
    if (M3ph_iw_nfft and M_iw_reg and density_matrix)
      chi3ph_iw_nfft = chi3_from_M3<Chan_t::PH>(M3ph_iw_nfft.value(), M_iw_reg.value(), G0_shift_iw_reg, density_matrix.value(), M_hartree.value());
    if (M3pp_iw_nfft_full and M_iw_reg and density_matrix)
      chi3pp_iw_nfft_full = chi3_from_M3<Chan_t::PP>(M3pp_iw_nfft_full.value(), M_iw_reg.value(), G0_shift_iw_reg, density_matrix.value(), M_hartree.value());
    if (M3ph_iw_nfft_full and M_iw_reg and density_matrix)
      chi3ph_iw_nfft_full = chi3_from_M3<Chan_t::PH>(M3ph_iw_nfft_full.value(), M_iw_reg.value(), G0_shift_iw_reg, density_matrix.value(), M_hartree.value());

    // Calculate chi2_iw from chi2_tau via DLR
    if (chi2pp_tau) chi2pp_iw = make_gf_dlr_imfreq(chi2pp_tau.value());
    if (chi2ph_tau) chi2ph_iw = make_gf_dlr_imfreq(chi2ph_tau.value());

    // Calculate chiAB_iw from chiAB_tau
    if (chiAB_tau) chiAB_iw = make_gf_dlr_imfreq(chiAB_tau.value());
  }

} // namespace triqs_ctint
