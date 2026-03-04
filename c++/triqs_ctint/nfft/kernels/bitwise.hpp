// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../common.hpp"

namespace triqs::utility::nfft {

  template <int Rank> struct kernel_bitwise_t {

    static constexpr int n_acc = 4; // ILP accumulators

    kernel_bitwise_t() = default;

    kernel_bitwise_t(shared_state_t<Rank> const &state) {
      unsigned long max_exponent = 0;
      for (int r = 0; r < Rank; ++r) {
        pow2_bits[r].resize(state.n_targets);
        for (int64_t d = 0; d < state.n_targets; ++d) {
          unsigned long exponent = odd_exponent_abs(state.target_n(r, d));
          max_exponent           = std::max(max_exponent, exponent);
          std::vector<int> bits;
          for (int k = 0; exponent > 0; ++k, exponent >>= 1)
            if (exponent & 1ul) bits.push_back(k);
          pow2_bits[r][d] = std::move(bits);
        }
      }
      num_pow2_levels = std::max(1, static_cast<int>(std::bit_width(max_exponent)));
      for (int r = 0; r < Rank; ++r) pow2_tbl[r].resize(num_pow2_levels, state.buf_size);
    }

    void execute(shared_state_t<Rank> &state) {
      double const pi_over_beta      = M_PI / state.beta;
      int64_t const buf_counter_simd = state.buf_counter & -simd_size;
      int64_t const stride           = state.buf_size;
      dcomplex *fiw_ptr              = state.fk_vec.data();

      // Phase 1: Build per-rank pow2 tables via sincos + repeated squaring
      poet::static_for<Rank>([&](const auto r) {
        dcomplex *tbl = pow2_tbl[r].data();

        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          using rbatch            = xsimd::batch<double>;
          auto [sin_vec, cos_vec] = xsimd::sincos(rbatch::load_unaligned(&state.x_arr(r, j)) * pi_over_beta);
          cbatch(cos_vec, sin_vec).store_unaligned(tbl + j);
        }
        for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
          double const theta = pi_over_beta * state.x_arr(r, j);
          tbl[j]             = dcomplex{std::cos(theta), std::sin(theta)};
        }

        for (int k = 1; k < num_pow2_levels; ++k) {
          dcomplex const *prev_row = tbl + (k - 1) * stride;
          dcomplex *cur_row        = tbl + k * stride;
          for (int j = 0; j < buf_counter_simd; j += simd_size) {
            cbatch prev = cbatch::load_unaligned(prev_row + j);
            (prev * prev).store_unaligned(cur_row + j);
          }
          for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
            dcomplex prev = prev_row[j];
            cur_row[j]    = prev * prev;
          }
        }
      });

      // Phase 2: accumulate with per-rank bit products
      std::array<dcomplex const *, Rank> tbl_base;
      poet::static_for<Rank>([&](const auto r) { tbl_base[r] = pow2_tbl[r].data(); });

      accumulate_targets_ilp<n_acc>(
         state.n_targets, buf_counter_simd, state.buf_counter, state.fx_arr.data(), fiw_ptr,
         [&](int64_t d, int j) -> cbatch {
           cbatch pow_prod;
           poet::static_for<Rank>([&](const auto r) {
             auto const *base = tbl_base[r];
             cbatch rank_pow(dcomplex{1.0, 0.0});
             for (int k : pow2_bits[r][d]) rank_pow *= cbatch::load_unaligned(base + k * stride + j);
             rank_pow = state.target_n(r, d) < 0 ? xsimd::conj(rank_pow) : rank_pow;
             pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
           });
           return pow_prod;
         },
         [&](int64_t d, int j) -> dcomplex {
           dcomplex pow_prod{1.0, 0.0};
           poet::static_for<Rank>([&](const auto r) {
             auto const *base = tbl_base[r];
             dcomplex rank_pow{1.0, 0.0};
             for (int k : pow2_bits[r][d]) rank_pow *= *(base + k * stride + j);
             rank_pow = state.target_n(r, d) < 0 ? std::conj(rank_pow) : rank_pow;
             pow_prod *= rank_pow;
           });
           return pow_prod;
         });
    }

    private:
    int num_pow2_levels = 0;
    std::array<nda::array<dcomplex, 2>, Rank> pow2_tbl;
    std::array<std::vector<std::vector<int>>, Rank> pow2_bits;
  };

} // namespace triqs::utility::nfft
