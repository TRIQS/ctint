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

  // Calibrate dispatch threshold for automatic mode by timing both NAF and Type3
  // at two buffer sizes and fitting a linear model to find the crossover.
  // Returns the buf_counter threshold below which NAF is preferred.
  template <int Rank>
  int calibrate_dispatch(shared_state_t<Rank> &state, kernel_naf_t<Rank> &naf_kernel, kernel_finufft_t<Rank> &finufft_kernel) {
    using clock = std::chrono::steady_clock;

    int const n_hi = std::min(4096, state.buf_size);
    int const n_lo = std::clamp(n_hi / 64, 1, 64);

    // Fill buffer with deterministic dummy data
    for (int j = 0; j < n_hi; ++j) {
      for (int r = 0; r < Rank; ++r) state.x_arr(r, j) = state.beta * static_cast<double>(j + 1) / (n_hi + 1);
      state.fx_arr[j] = dcomplex(1.0, 0.0);
    }

    auto measure_naf = [&](int n) {
      state.buf_counter = n;
      state.fk_vec      = 0;
      auto t0           = clock::now();
      naf_kernel.execute(state);
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
    measure_naf(n_lo);

    // Measure at two buffer sizes (take best of 3 reps)
    constexpr int n_reps = 3;
    auto best_of = [&](auto &&fn, int n) {
      double best = std::numeric_limits<double>::max();
      for (int rep = 0; rep < n_reps; ++rep) best = std::min(best, fn(n));
      return best;
    };

    double naf_lo = best_of(measure_naf, n_lo);
    double naf_hi = best_of(measure_naf, n_hi);
    double t3_lo  = best_of(measure_type3, n_lo);
    double t3_hi  = best_of(measure_type3, n_hi);

    // Linear fit: time(B) = intercept + slope * B
    double naf_slope     = (naf_hi - naf_lo) / (n_hi - n_lo);
    double naf_intercept = naf_lo - naf_slope * n_lo;
    double t3_slope      = (t3_hi - t3_lo) / (n_hi - n_lo);
    double t3_intercept  = t3_lo - t3_slope * n_lo;

    int threshold;
    // Crossover: naf_intercept + naf_slope * B = t3_intercept + t3_slope * B
    if (naf_slope > t3_slope && t3_intercept > naf_intercept)
      threshold = static_cast<int>((t3_intercept - naf_intercept) / (naf_slope - t3_slope));
    else if (naf_slope <= t3_slope)
      threshold = state.buf_size + 1; // NAF always wins
    else
      threshold = 0; // Type3 always wins

    threshold = std::clamp(threshold, 0, state.buf_size + 1);

    // Clean up
    state.fk_vec      = 0;
    state.buf_counter = 0;

    return threshold;
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

    double dir_slope     = (dir_hi - dir_lo) / (n_hi - n_lo);
    double dir_intercept = dir_lo - dir_slope * n_lo;
    double t1_slope      = (t1_hi - t1_lo) / (n_hi - n_lo);
    double t1_intercept  = t1_lo - t1_slope * n_lo;

    int threshold;
    if (dir_slope > t1_slope && t1_intercept > dir_intercept)
      threshold = static_cast<int>((t1_intercept - dir_intercept) / (dir_slope - t1_slope));
    else if (dir_slope <= t1_slope)
      threshold = state.buf_size + 1; // direct always wins
    else
      threshold = 0; // FINUFFT type1 always wins

    threshold = std::clamp(threshold, 0, state.buf_size + 1);

    // Clean up
    restore_state();
    state.fk_vec      = 0;
    state.buf_counter = 0;

    return threshold;
  }

} // namespace triqs::utility::nfft
