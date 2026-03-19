// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "common.hpp"
#include "kernels/direct_type1.hpp"
#include "kernels/finufft.hpp"
#include "kernels/naf.hpp"
#include <chrono>

namespace triqs::utility::nfft {

  // ---- Calibration helpers ----

  // Linear-fit crossover: given timings of two kernels at n_lo and n_hi,
  // find the buffer size where they cross. Returns threshold clamped to [0, buf_size+1].
  inline int linear_crossover(double a_lo, double a_hi, double b_lo, double b_hi, int n_lo, int n_hi, int buf_size) {
    double a_slope     = (a_hi - a_lo) / (n_hi - n_lo);
    double a_intercept = a_lo - a_slope * n_lo;
    double b_slope     = (b_hi - b_lo) / (n_hi - n_lo);
    double b_intercept = b_lo - b_slope * n_lo;

    int threshold;
    if (a_slope > b_slope && b_intercept > a_intercept)
      threshold = static_cast<int>((b_intercept - a_intercept) / (a_slope - b_slope));
    else if (a_slope <= b_slope)
      threshold = buf_size + 1; // a always wins
    else
      threshold = 0; // b always wins

    return std::clamp(threshold, 0, buf_size + 1);
  }

  // Fill buffer with deterministic dummy tau data for calibration
  template <int Rank> void fill_dummy_tau(shared_state_t<Rank> &state, int n) {
    for (int j = 0; j < n; ++j) {
      for (int r = 0; r < Rank; ++r) state.x_arr(r, j) = state.beta * static_cast<double>(j + 1) / (n + 1);
      state.fx_arr[j] = dcomplex(1.0, 0.0);
    }
  }

  // Minimum of n_reps evaluations of fn(n)
  template <typename F> double best_of_n(F &&fn, int n, int n_reps = 3) {
    double best = std::numeric_limits<double>::max();
    for (int rep = 0; rep < n_reps; ++rep) best = std::min(best, fn(n));
    return best;
  }

  // Save/restore calibration state (FINUFFT transforms modify x_arr/fx_arr in-place)
  template <int Rank> struct state_saver_t {
    shared_state_t<Rank> &state;
    int n;
    nda::array<double, 2> x_save;
    nda::vector<dcomplex> fx_save;

    state_saver_t(shared_state_t<Rank> &state_, int n_) : state(state_), n(n_), x_save(Rank, n_), fx_save(n_) { save(); }
    ~state_saver_t() { cleanup(); }

    void save() {
      x_save  = state.x_arr(nda::range(0, Rank), nda::range(0, n));
      fx_save = state.fx_arr(nda::range(0, n));
    }

    void restore() {
      state.x_arr(nda::range(0, Rank), nda::range(0, n)) = x_save;
      state.fx_arr(nda::range(0, n))                     = fx_save;
    }

    void cleanup() {
      restore();
      state.fk_vec      = 0;
      state.buf_counter = 0;
    }
  };

  // ---- Calibration functions ----

  // Selection plan for automatic non-uniform dispatch.
  // There are two independent decisions:
  // 1. Direct vs FINUFFT: below direct_vs_finufft_threshold use a direct kernel, above it use FINUFFT.
  // 2. Inside the FINUFFT region: switch between type3 and type1_gather at finufft_path_switch_threshold.
  struct nonuniform_dispatch_plan_t {
    int direct_vs_finufft_threshold;
    int finufft_path_switch_threshold;
    bool select_direct_type1;
    bool use_type3_below_switch_threshold;
    double direct_path_lo_time;
    double direct_path_hi_time;
  };

  // Calibrate dispatch threshold for automatic mode with non-uniform targets.
  // Compares direct kernels (direct_type1 vs sparse direct kernel) against both
  // FINUFFT paths (type1_gather vs type3), then picks the faster direct kernel
  // and the faster FINUFFT path, and computes their crossover.
  template <int TolDigits = 12, int Rank, typename SparseDirectKernel, typename FinufftKernel>
  nonuniform_dispatch_plan_t calibrate_dispatch_nonuniform(shared_state_t<Rank> &state, kernel_direct_type1_t<Rank> &direct_type1_kernel,
                                                           SparseDirectKernel &sparse_direct_kernel, FinufftKernel &finufft_kernel) {
    using clock = std::chrono::steady_clock;

    int const n_hi = std::min(4096, state.buf_size);
    // Rank-1 auto workloads frequently live in the n=8..128 regime, so calibrating only
    // from n=64 upward can switch to the FINUFFT side too early around the crossover.
    int const n_lo = [] {
      if constexpr (Rank == 1)
        return 8;
      else
        return 64;
    }();

    fill_dummy_tau(state, n_hi);
    state_saver_t saver(state, n_hi);

    auto measure_kernel = [&](auto &kernel, int n) {
      saver.restore();
      state.buf_counter = n;
      state.fk_vec      = 0;
      auto t0           = clock::now();
      kernel.execute(state);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    nda::array<dcomplex, 1> dummy_fiw(state.n_targets);

    auto measure_type1_gather = [&](int n) {
      saver.restore();
      state.buf_counter = n;
      apply_type1_coord_transform<TolDigits>(state, n);
      dummy_fiw = 0;
      auto t0   = clock::now();
      finufft_kernel.execute_type1_gather(state, dummy_fiw);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    auto measure_type3 = [&](int n) {
      saver.restore();
      state.buf_counter = n;
      dummy_fiw         = 0;
      auto t0           = clock::now();
      finufft_kernel.execute_type3(state, dummy_fiw);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    // Warmup
    measure_type1_gather(n_lo);
    measure_type3(n_lo);
    measure_kernel(sparse_direct_kernel, n_lo);
    measure_kernel(direct_type1_kernel, n_lo);

    // Pick faster direct kernel at n_lo, measure winner at n_hi
    double dt1_lo    = best_of_n([&](int n) { return measure_kernel(direct_type1_kernel, n); }, n_lo);
    double sparse_lo = best_of_n([&](int n) { return measure_kernel(sparse_direct_kernel, n); }, n_lo);
    bool const select_direct_type1 = dt1_lo <= sparse_lo;
    double dir_lo                  = select_direct_type1 ? dt1_lo : sparse_lo;
    double dir_hi                  = select_direct_type1 ? best_of_n([&](int n) { return measure_kernel(direct_type1_kernel, n); }, n_hi) :
                                                           best_of_n([&](int n) { return measure_kernel(sparse_direct_kernel, n); }, n_hi);

    // Measure both FINUFFT paths. Automatic mode can dispatch between them, so calibrate
    // both the per-n FINUFFT crossover and the direct-vs-best-FINUFFT crossover.
    double tg_lo = best_of_n(measure_type1_gather, n_lo);
    double tg_hi = best_of_n(measure_type1_gather, n_hi);
    double t3_lo = best_of_n(measure_type3, n_lo);
    double t3_hi = best_of_n(measure_type3, n_hi);

    bool const use_type3_below_switch_threshold = t3_lo <= tg_lo;
    int const finufft_path_switch_threshold =
       use_type3_below_switch_threshold ? linear_crossover(t3_lo, t3_hi, tg_lo, tg_hi, n_lo, n_hi, state.buf_size) :
                                          linear_crossover(tg_lo, tg_hi, t3_lo, t3_hi, n_lo, n_hi, state.buf_size);

    double finufft_lo = std::min(t3_lo, tg_lo);
    double finufft_hi = std::min(t3_hi, tg_hi);

    int direct_vs_finufft_threshold = linear_crossover(dir_lo, dir_hi, finufft_lo, finufft_hi, n_lo, n_hi, state.buf_size);
    if constexpr (Rank == 1) direct_vs_finufft_threshold = std::max(direct_vs_finufft_threshold, 96);

    saver.cleanup();
    return {.direct_vs_finufft_threshold = direct_vs_finufft_threshold,
            .finufft_path_switch_threshold = finufft_path_switch_threshold,
            .select_direct_type1 = select_direct_type1,
            .use_type3_below_switch_threshold = use_type3_below_switch_threshold,
            .direct_path_lo_time = dir_lo,
            .direct_path_hi_time = dir_hi};
  }

  // Calibrate dispatch threshold for type1 automatic mode (direct_type1 vs FINUFFT type1).
  // The direct path uses raw tau, while the FINUFFT path needs coordinate transformation.
  template <int TolDigits = 12, int Rank, typename FinufftKernel>
  int calibrate_dispatch_type1(shared_state_t<Rank> &state, kernel_direct_type1_t<Rank> &direct_kernel, FinufftKernel &finufft_kernel,
                               nda::array<dcomplex, Rank> &fk_arr, int common_factor) {
    using clock = std::chrono::steady_clock;

    int const n_hi = std::min(4096, state.buf_size);
    int const n_lo = std::clamp(n_hi / 64, 1, 64);

    fill_dummy_tau(state, n_hi);
    state_saver_t saver(state, n_hi);

    auto measure_direct = [&](int n) {
      saver.restore();
      state.buf_counter = n;
      state.fk_vec      = 0;
      auto t0           = clock::now();
      direct_kernel.execute(state);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    nda::array<dcomplex, Rank> dummy_fiw(fk_arr.shape());

    auto measure_type1 = [&](int n) {
      saver.restore();
      state.buf_counter = n;
      apply_type1_coord_transform<TolDigits>(state, n);
      dummy_fiw = 0;
      auto t0   = clock::now();
      finufft_kernel.execute_type1(state, dummy_fiw, fk_arr, common_factor);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    // Warmup
    measure_type1(n_lo);
    measure_direct(n_lo);

    double dir_lo = best_of_n(measure_direct, n_lo);
    double dir_hi = best_of_n(measure_direct, n_hi);
    double t1_lo  = best_of_n(measure_type1, n_lo);
    double t1_hi  = best_of_n(measure_type1, n_hi);

    int threshold = linear_crossover(dir_lo, dir_hi, t1_lo, t1_hi, n_lo, n_hi, state.buf_size);

    saver.cleanup();
    return threshold;
  }

} // namespace triqs::utility::nfft
