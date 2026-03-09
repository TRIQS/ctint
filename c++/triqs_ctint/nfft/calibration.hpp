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

  // Calibrate dispatch threshold for automatic mode (non-uniform targets).
  // Compares direct_type1, NAF, and FINUFFT type3 kernels:
  // 1. Pick the faster direct kernel (direct_type1 vs NAF) at a representative buffer size
  // 2. Calibrate the winner against FINUFFT type3 using linear fit
  // Returns {threshold, use_direct_type1}: threshold below which the winning direct kernel is preferred.
  template <int Rank>
  std::pair<int, bool> calibrate_dispatch(shared_state_t<Rank> &state, kernel_direct_type1_t<Rank> &direct_type1_kernel,
                                          kernel_naf_t<Rank> &naf_kernel, kernel_finufft_t<Rank> &finufft_kernel) {
    using clock = std::chrono::steady_clock;

    int const n_hi = std::min(4096, state.buf_size);
    int const n_lo = std::clamp(n_hi / 64, 1, 64);

    // Fill buffer with deterministic dummy data
    for (int j = 0; j < n_hi; ++j) {
      for (int r = 0; r < Rank; ++r) state.x_arr(r, j) = state.beta * static_cast<double>(j + 1) / (n_hi + 1);
      state.fx_arr[j] = dcomplex(1.0, 0.0);
    }

    auto measure_kernel = [&](auto &kernel, int n) {
      state.buf_counter = n;
      state.fk_vec      = 0;
      auto t0           = clock::now();
      kernel.execute(state);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    auto measure_type3 = [&](int n) {
      state.buf_counter = n;
      auto t0           = clock::now();
      finufft_kernel.set_pts_type3(state);
      check_finufft(finufft_execute(finufft_kernel.get_plan(), state.fx_arr.data(), state.fk_vec.data()));
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    // Warmup to avoid cold-start bias
    measure_type3(n_lo);
    measure_kernel(naf_kernel, n_lo);
    measure_kernel(direct_type1_kernel, n_lo);

    // Measure at two buffer sizes (take best of 3 reps)
    constexpr int n_reps = 3;
    auto best_of = [&](auto &&fn, int n) {
      double best = std::numeric_limits<double>::max();
      for (int rep = 0; rep < n_reps; ++rep) best = std::min(best, fn(n));
      return best;
    };

    // Pick the faster direct kernel at n_lo, then only measure winner at n_hi.
    // Note: winner is chosen at n_lo only; both kernels are O(n * n_targets) so
    // relative ranking is expected to hold across buffer sizes.
    double dt1_lo = best_of([&](int n) { return measure_kernel(direct_type1_kernel, n); }, n_lo);
    double naf_lo = best_of([&](int n) { return measure_kernel(naf_kernel, n); }, n_lo);
    bool use_dt1  = dt1_lo <= naf_lo;
    double dir_lo = use_dt1 ? dt1_lo : naf_lo;
    double dir_hi = use_dt1 ? best_of([&](int n) { return measure_kernel(direct_type1_kernel, n); }, n_hi)
                            : best_of([&](int n) { return measure_kernel(naf_kernel, n); }, n_hi);
    double t3_lo  = best_of(measure_type3, n_lo);
    double t3_hi  = best_of(measure_type3, n_hi);

    int threshold = linear_crossover(dir_lo, dir_hi, t3_lo, t3_hi, n_lo, n_hi, state.buf_size);

    // Clean up
    state.fk_vec      = 0;
    state.buf_counter = 0;

    return {threshold, use_dt1};
  }

  // Calibrate dispatch threshold for type1 automatic mode (direct_type1 vs FINUFFT type1).
  // The direct path uses raw tau, while the FINUFFT path needs coordinate transformation.
  // We time the full alternative: direct_type1 execute vs (coord transform + FINUFFT type1 execute).
  template <int Rank>
  int calibrate_dispatch_type1(shared_state_t<Rank> &state, kernel_direct_type1_t<Rank> &direct_kernel, kernel_finufft_t<Rank> &finufft_kernel,
                               nda::array<dcomplex, Rank> &fk_arr, int common_factor) {
    using clock = std::chrono::steady_clock;

    int const n_hi = std::min(4096, state.buf_size);
    int const n_lo = std::clamp(n_hi / 64, 1, 64);

    // Fill buffer with deterministic dummy raw tau data
    for (int j = 0; j < n_hi; ++j) {
      for (int r = 0; r < Rank; ++r) state.x_arr(r, j) = state.beta * static_cast<double>(j + 1) / (n_hi + 1);
      state.fx_arr[j] = dcomplex(1.0, 0.0);
    }

    // Save raw tau for reuse (FINUFFT transform modifies x_arr/fx_arr in-place)
    nda::array<double, 2> x_arr_save(state.x_arr.shape());
    nda::vector<dcomplex> fx_arr_save(state.fx_arr.shape()[0]);

    auto save_state = [&]() {
      x_arr_save(nda::range(0, Rank), nda::range(0, n_hi)) = state.x_arr(nda::range(0, Rank), nda::range(0, n_hi));
      for (int j = 0; j < n_hi; ++j) fx_arr_save[j] = state.fx_arr[j];
    };

    auto restore_state = [&]() {
      state.x_arr(nda::range(0, Rank), nda::range(0, n_hi)) = x_arr_save(nda::range(0, Rank), nda::range(0, n_hi));
      for (int j = 0; j < n_hi; ++j) state.fx_arr[j] = fx_arr_save[j];
    };

    save_state();

    auto measure_direct = [&](int n) {
      restore_state();
      state.buf_counter = n;
      state.fk_vec      = 0;
      auto t0           = clock::now();
      direct_kernel.execute(state);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    // Dummy array view for execute_type1 output (we don't care about correctness, only timing)
    nda::array<dcomplex, Rank> dummy_fiw(fk_arr.shape());

    auto measure_type1 = [&](int n) {
      restore_state();
      state.buf_counter = n;
      // Transform coordinates (part of the FINUFFT path cost)
      double const inv_beta = 1.0 / state.beta;
      for (int j = 0; j < n; ++j) {
        double tau_sum = 0.0;
        for (int r = 0; r < Rank; ++r) {
          double tau = state.x_arr(r, j);
          tau_sum += tau;
          state.x_arr(r, j) = 2 * M_PI * (tau * inv_beta - 0.5);
        }
        state.fx_arr[j] *= std::exp(dcomplex(0, M_PI * tau_sum * inv_beta));
      }
      dummy_fiw = 0;
      auto t0 = clock::now();
      finufft_kernel.execute_type1(state, dummy_fiw, fk_arr, common_factor);
      return std::chrono::duration<double>(clock::now() - t0).count();
    };

    // Warmup
    measure_type1(n_lo);
    measure_direct(n_lo);

    constexpr int n_reps = 3;
    auto best_of = [&](auto &&fn, int n) {
      double best = std::numeric_limits<double>::max();
      for (int rep = 0; rep < n_reps; ++rep) best = std::min(best, fn(n));
      return best;
    };

    double dir_lo = best_of(measure_direct, n_lo);
    double dir_hi = best_of(measure_direct, n_hi);
    double t1_lo  = best_of(measure_type1, n_lo);
    double t1_hi  = best_of(measure_type1, n_hi);

    int threshold = linear_crossover(dir_lo, dir_hi, t1_lo, t1_hi, n_lo, n_hi, state.buf_size);

    // Clean up
    restore_state();
    state.fk_vec      = 0;
    state.buf_counter = 0;

    return threshold;
  }

} // namespace triqs::utility::nfft
