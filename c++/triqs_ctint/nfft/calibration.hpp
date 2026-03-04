// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "common.hpp"
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

} // namespace triqs::utility::nfft
