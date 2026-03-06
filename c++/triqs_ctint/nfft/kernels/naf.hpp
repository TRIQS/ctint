// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../common.hpp"
#include <unordered_map>

namespace triqs::utility::nfft {

  template <int Rank> struct kernel_naf_t {

    kernel_naf_t() = default;

    kernel_naf_t(shared_state_t<Rank> const &state, int buf_size) {
      unsigned long max_exponent = 0;
      target_map.resize(state.n_targets * Rank);
      for (int r = 0; r < Rank; ++r) {
        std::unordered_map<unsigned long, int> exp_to_idx;
        digit_offsets[r].push_back(0);
        for (int64_t d = 0; d < state.n_targets; ++d) {
          unsigned long exp   = odd_exponent_abs(state.target_n(r, d));
          max_exponent        = std::max(max_exponent, exp);
          auto [it, inserted] = exp_to_idx.try_emplace(exp, static_cast<int>(exp_to_idx.size()));
          if (inserted) {
            auto digits = compute_naf(exp);
            digits_flat[r].insert(digits_flat[r].end(), digits.begin(), digits.end());
            digit_offsets[r].push_back(static_cast<int>(digits_flat[r].size()));
          }
          bool needs_conj          = (state.target_n(r, d) < 0);
          target_map[d * Rank + r] = 2 * it->second + (needs_conj ? 1 : 0);
        }
        n_unique[r] = static_cast<int>(exp_to_idx.size());
      }
      num_pow2_levels = std::max(1, static_cast<int>(std::bit_width(max_exponent)) + 1);
      for (int r = 0; r < Rank; ++r) {
        pow2_tbl[r].resize(num_pow2_levels, buf_size);
        uq_simd_buf[r].resize(2 * n_unique[r]);
        uq_scalar_buf[r].resize(2 * n_unique[r]);
      }
      sums_buf.resize(state.n_targets);
    }

    void execute(shared_state_t<Rank> &state) {
      double const pi_over_beta      = M_PI / state.beta;
      int64_t const buf_counter_simd = state.buf_counter & -simd_size;
      int64_t const stride           = state.buf_size;

      // Phase 1: Build per-rank pow2 tables via repeated squaring
      poet::static_for<Rank>([&](const auto r) {
        dcomplex *tbl = pow2_tbl[r].data();

        // Level 0: z_r = exp(i*pi*tau_r/beta)
        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          using rbatch            = xsimd::batch<double>;
          auto [sin_vec, cos_vec] = xsimd::sincos(rbatch::load_unaligned(&state.x_arr(r, j)) * pi_over_beta);
          cbatch(cos_vec, sin_vec).store_unaligned(tbl + j);
        }
        for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
          double const theta = pi_over_beta * state.x_arr(r, j);
          tbl[j]             = dcomplex{std::cos(theta), std::sin(theta)};
        }

        // Higher levels: squaring
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

      // Phase 2: Factored unique-exponent accumulation.
      execute_phase2(state);
    }

    private:
    int num_pow2_levels = 0;
    std::array<nda::array<dcomplex, 2>, Rank> pow2_tbl;
    std::array<int, Rank> n_unique{};
    std::array<std::vector<int>, Rank> digits_flat;
    std::array<std::vector<int>, Rank> digit_offsets;
    std::vector<int> target_map;
    mutable std::array<std::vector<xsimd::batch<dcomplex>>, Rank> uq_simd_buf;
    mutable std::array<std::vector<dcomplex>, Rank> uq_scalar_buf;
    mutable std::vector<xsimd::batch<dcomplex>> sums_buf;

    void execute_phase2(shared_state_t<Rank> &state) {
      int64_t const buf_counter_simd = state.buf_counter & -simd_size;
      int64_t const stride           = state.buf_size;
      dcomplex *fiw_ptr              = state.fk_vec.data();

      std::array<dcomplex const *, Rank> tbl_base;
      poet::static_for<Rank>([&](const auto r) { tbl_base[r] = pow2_tbl[r].data(); });

      int const *map_ptr                = target_map.data();
      constexpr int map_stride          = Rank;
      constexpr int64_t l2_bytes        = 2 * 1024 * 1024;
      int64_t const table_bytes_per_src = Rank * num_pow2_levels * static_cast<int64_t>(sizeof(dcomplex));
      bool const use_blocking           = buf_counter_simd * table_bytes_per_src > l2_bytes;

      auto compute_unique_simd = [&](int r, int u, int j) -> cbatch {
        int const *digs = digits_flat[r].data() + digit_offsets[r][u];
        int n_digits    = digit_offsets[r][u + 1] - digit_offsets[r][u];
        auto const *base = tbl_base[r];

        int d0          = digs[0];
        int row0        = d0 >= 0 ? d0 : -(d0 + 1);
        cbatch rank_pow = cbatch::load_unaligned(base + row0 * stride + j);
        if (d0 < 0) rank_pow = xsimd::conj(rank_pow);

        for (int i = 1; i < n_digits; ++i) {
          int di     = digs[i];
          int row    = di >= 0 ? di : -(di + 1);
          cbatch val = cbatch::load_unaligned(base + row * stride + j);
          rank_pow *= di >= 0 ? val : xsimd::conj(val);
        }
        return rank_pow;
      };

      auto compute_unique_scalar = [&](int r, int u, int j) -> dcomplex {
        int const *digs  = digits_flat[r].data() + digit_offsets[r][u];
        int n_digits     = digit_offsets[r][u + 1] - digit_offsets[r][u];
        auto const *base = tbl_base[r];

        int d0            = digs[0];
        int row0          = d0 >= 0 ? d0 : -(d0 + 1);
        dcomplex rank_pow = *(base + row0 * stride + j);
        if (d0 < 0) rank_pow = std::conj(rank_pow);

        for (int i = 1; i < n_digits; ++i) {
          int di       = digs[i];
          int row      = di >= 0 ? di : -(di + 1);
          dcomplex val = *(base + row * stride + j);
          rank_pow *= di >= 0 ? val : std::conj(val);
        }
        return rank_pow;
      };

      auto accumulate_simd_range = [&](int j_begin, int j_end) {
        // Restrict-qualified local pointers break aliasing assumptions that prevent
        // the compiler from hoisting data pointers out of the inner target loop.
        cbatch *__restrict__ sp            = sums_buf.data();
        cbatch *__restrict__ uq0           = uq_simd_buf[0].data();
        int const *__restrict__ mp         = map_ptr;
        int64_t const n_tgt                = state.n_targets;
        dcomplex const *__restrict__ fxp   = state.fx_arr.data();

        [[maybe_unused]] cbatch *__restrict__ uq1 = nullptr;
        if constexpr (Rank >= 2) uq1 = uq_simd_buf[1].data();

        auto fill_unique = [&](cbatch *__restrict__ uq, auto r, int j_idx) {
          for (int u = 0; u < n_unique[r]; ++u) {
            cbatch val     = compute_unique_simd(r, u, j_idx);
            uq[2 * u]     = val;
            uq[2 * u + 1] = xsimd::conj(val);
          }
        };

        for (int j = j_begin; j < j_end; j += simd_size) {
          fill_unique(uq0, std::integral_constant<int, 0>{}, j);
          if constexpr (Rank >= 2) fill_unique(uq1, std::integral_constant<int, 1>{}, j);
          if constexpr (Rank > 2) {
            for (int r = 2; r < Rank; ++r) fill_unique(uq_simd_buf[r].data(), r, j);
          }
          cbatch fj = cbatch::load_unaligned(fxp + j);
          for (int64_t d = 0; d < n_tgt; ++d) {
            int const *info = mp + d * map_stride;
            cbatch pow      = uq0[info[0]];
            if constexpr (Rank >= 2) pow *= uq1[info[1]];
            if constexpr (Rank > 2) {
              for (int r = 2; r < Rank; ++r) pow *= uq_simd_buf[r][info[r]];
            }
            sp[d] = xsimd::fma(fj, pow, sp[d]);
          }
        }
      };

      std::fill(sums_buf.begin(), sums_buf.end(), cbatch(dcomplex{0, 0}));
      if (use_blocking) {
        constexpr int source_block = 128;
        for (int jb = 0; jb < buf_counter_simd; jb += source_block)
          accumulate_simd_range(jb, std::min(jb + source_block, static_cast<int>(buf_counter_simd)));
      } else {
        accumulate_simd_range(0, static_cast<int>(buf_counter_simd));
      }
      for (int64_t d = 0; d < state.n_targets; ++d) fiw_ptr[d] += xsimd::reduce_add(sums_buf[d]);

      // Scalar tail
      for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
        poet::static_for<Rank>([&](const auto r) {
          for (int u = 0; u < n_unique[r]; ++u) {
            dcomplex val              = compute_unique_scalar(r, u, j);
            uq_scalar_buf[r][2 * u]     = val;
            uq_scalar_buf[r][2 * u + 1] = std::conj(val);
          }
        });
        dcomplex fj = state.fx_arr[j];
        for (int64_t d = 0; d < state.n_targets; ++d) {
          int const *info = map_ptr + d * map_stride;
          dcomplex pow    = uq_scalar_buf[0][info[0]];
          for (int r = 1; r < Rank; ++r) pow *= uq_scalar_buf[r][info[r]];
          fiw_ptr[d] += fj * pow;
        }
      }
    }
  };

} // namespace triqs::utility::nfft
