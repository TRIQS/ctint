// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../common.hpp"
#include <unordered_map>

namespace triqs::utility::nfft {

  template <int Rank> struct kernel_naf_t {

    static constexpr int n_acc = 4; // ILP accumulators

    kernel_naf_t() = default;

    kernel_naf_t(shared_state_t<Rank> const &state, int buf_size) {
      unsigned long max_exponent = 0;
      for (int r = 0; r < Rank; ++r) {
        digit_offsets[r].resize(state.n_targets + 1);
        digit_offsets[r][0] = 0;
        for (int64_t d = 0; d < state.n_targets; ++d) {
          unsigned long exp = odd_exponent_abs(state.target_n(r, d));
          max_exponent      = std::max(max_exponent, exp);
          auto digits       = compute_naf(exp);
          digits_flat[r].insert(digits_flat[r].end(), digits.begin(), digits.end());
          digit_offsets[r][d + 1] = static_cast<int>(digits_flat[r].size());
        }
      }
      num_pow2_levels = std::max(1, static_cast<int>(std::bit_width(max_exponent)) + 1);
      for (int r = 0; r < Rank; ++r) pow2_tbl[r].resize(num_pow2_levels, buf_size);

      // For Rank >= 2, build factored unique-exponent data.
      if constexpr (Rank >= 2) {
        int total_unique = 0;
        target_map.resize(state.n_targets * (Rank + 1));
        for (int r = 0; r < Rank; ++r) {
          std::unordered_map<unsigned long, int> exp_to_idx;
          uniq_digit_offsets[r].push_back(0);
          for (int64_t d = 0; d < state.n_targets; ++d) {
            unsigned long exp   = odd_exponent_abs(state.target_n(r, d));
            auto [it, inserted] = exp_to_idx.try_emplace(exp, static_cast<int>(exp_to_idx.size()));
            if (inserted) {
              auto digits = compute_naf(exp);
              uniq_digits_flat[r].insert(uniq_digits_flat[r].end(), digits.begin(), digits.end());
              uniq_digit_offsets[r].push_back(static_cast<int>(uniq_digits_flat[r].size()));
            }
            target_map[d * (Rank + 1) + r] = it->second;
          }
          n_unique[r] = static_cast<int>(exp_to_idx.size());
          total_unique += n_unique[r];
        }
        // Store conjugation bitmask per target (bit r set if target_n(r,d) < 0)
        for (int64_t d = 0; d < state.n_targets; ++d) {
          int conj_mask = 0;
          for (int r = 0; r < Rank; ++r)
            if (state.target_n(r, d) < 0) conj_mask |= (1 << r);
          target_map[d * (Rank + 1) + Rank] = conj_mask;
        }
        use_factored = total_unique < state.n_targets * Rank;
        if (use_factored) {
          for (int r = 0; r < Rank; ++r) {
            uq_simd_buf[r].resize(n_unique[r]);
            uq_scalar_buf[r].resize(n_unique[r]);
          }
          sums_buf.resize(state.n_targets);
        }
      }
    }

    void execute(shared_state_t<Rank> &state) {
      double const pi_over_beta      = M_PI / state.beta;
      int64_t const buf_counter_simd = state.buf_counter & -simd_size;
      int64_t const stride           = state.buf_size;
      dcomplex *fiw_ptr              = state.fk_vec.data();

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

      // Phase 2: Dispatch to factored path for Rank >= 2 when beneficial.
      if constexpr (Rank >= 2) {
        if (use_factored) {
          execute_phase2_factored(state);
          return;
        }
      }

      // Phase 2 (non-factored): Source-blocked accumulation for cache locality.
      std::array<dcomplex const *, Rank> tbl_base;
      poet::static_for<Rank>([&](const auto r) { tbl_base[r] = pow2_tbl[r].data(); });

      auto compute_simd_pow = [&](int64_t d, int j) -> cbatch {
        cbatch pow_prod;
        poet::static_for<Rank>([&](const auto r) {
          int const *digs  = digits_flat[r].data() + digit_offsets[r][d];
          int n_digits     = digit_offsets[r][d + 1] - digit_offsets[r][d];
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

          rank_pow = (state.target_n(r, d) < 0) ? xsimd::conj(rank_pow) : rank_pow;
          pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
        });
        return pow_prod;
      };

      auto compute_scalar_pow = [&](int64_t d, int j) -> dcomplex {
        dcomplex pow_prod{1.0, 0.0};
        poet::static_for<Rank>([&](const auto r) {
          int const *digs  = digits_flat[r].data() + digit_offsets[r][d];
          int n_digits     = digit_offsets[r][d + 1] - digit_offsets[r][d];
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

          rank_pow = (state.target_n(r, d) < 0) ? std::conj(rank_pow) : rank_pow;
          pow_prod *= rank_pow;
        });
        return pow_prod;
      };

      // Use source blocking when the table exceeds L2 cache, otherwise use unblocked ILP.
      constexpr int64_t l2_bytes           = 2 * 1024 * 1024;
      int64_t const table_bytes_per_source = Rank * num_pow2_levels * static_cast<int64_t>(sizeof(dcomplex));
      bool const use_blocking              = buf_counter_simd * table_bytes_per_source > l2_bytes;

      if (use_blocking) {
        constexpr int source_block   = 128;
        int64_t const n_targets_main = (state.n_targets / n_acc) * n_acc;

        for (int jb = 0; jb < buf_counter_simd; jb += source_block) {
          int const j_end = std::min(jb + source_block, static_cast<int>(buf_counter_simd));

          int64_t d = 0;
          for (; d < n_targets_main; d += n_acc) {
            std::array<cbatch, n_acc> local_sums;
            poet::static_for<n_acc>([&](const auto i) { local_sums[i] = cbatch(dcomplex{0, 0}); });

            for (int j = jb; j < j_end; j += simd_size) {
              cbatch fj = cbatch::load_unaligned(state.fx_arr.data() + j);
              poet::static_for<n_acc>(
                 [&](const auto i) { local_sums[i] = xsimd::fma(fj, compute_simd_pow(d + i, j), local_sums[i]); });
            }

            poet::static_for<n_acc>([&](const auto i) { fiw_ptr[d + i] += xsimd::reduce_add(local_sums[i]); });
          }
          for (; d < state.n_targets; ++d) {
            cbatch local_sum(dcomplex{0, 0});
            for (int j = jb; j < j_end; j += simd_size)
              local_sum = xsimd::fma(cbatch::load_unaligned(state.fx_arr.data() + j), compute_simd_pow(d, j), local_sum);
            fiw_ptr[d] += xsimd::reduce_add(local_sum);
          }
        }

        // Scalar tail
        for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
          dcomplex fj = state.fx_arr[j];
          for (int64_t d = 0; d < state.n_targets; ++d) fiw_ptr[d] += fj * compute_scalar_pow(d, j);
        }
      } else {
        accumulate_targets_ilp<n_acc>(state.n_targets, buf_counter_simd, state.buf_counter, state.fx_arr.data(), fiw_ptr,
                                      compute_simd_pow, compute_scalar_pow);
      }
    }

    private:
    int num_pow2_levels = 0;
    std::array<nda::array<dcomplex, 2>, Rank> pow2_tbl;
    std::array<std::vector<int>, Rank> digits_flat;
    std::array<std::vector<int>, Rank> digit_offsets;

    // Rank >= 2 factored NAF data
    bool use_factored = false;
    std::array<int, Rank> n_unique{};
    std::array<std::vector<int>, Rank> uniq_digits_flat;
    std::array<std::vector<int>, Rank> uniq_digit_offsets;
    std::vector<int> target_map;
    mutable std::array<std::vector<xsimd::batch<dcomplex>>, Rank> uq_simd_buf;
    mutable std::array<std::vector<dcomplex>, Rank> uq_scalar_buf;
    mutable std::vector<xsimd::batch<dcomplex>> sums_buf;

    void execute_phase2_factored(shared_state_t<Rank> &state)
      requires(Rank >= 2)
    {
      int64_t const buf_counter_simd = state.buf_counter & -simd_size;
      int64_t const stride           = state.buf_size;
      dcomplex *fiw_ptr              = state.fk_vec.data();

      std::array<dcomplex const *, Rank> tbl_base;
      poet::static_for<Rank>([&](const auto r) { tbl_base[r] = pow2_tbl[r].data(); });

      int const *map_ptr                = target_map.data();
      constexpr int map_stride          = Rank + 1;
      constexpr int64_t l2_bytes        = 2 * 1024 * 1024;
      int64_t const table_bytes_per_src = Rank * num_pow2_levels * static_cast<int64_t>(sizeof(dcomplex));
      bool const use_blocking           = buf_counter_simd * table_bytes_per_src > l2_bytes;

      auto compute_unique_simd = [&](int r, int u, int j) -> cbatch {
        int const *digs = uniq_digits_flat[r].data() + uniq_digit_offsets[r][u];
        int n_digits    = uniq_digit_offsets[r][u + 1] - uniq_digit_offsets[r][u];
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
        int const *digs  = uniq_digits_flat[r].data() + uniq_digit_offsets[r][u];
        int n_digits     = uniq_digit_offsets[r][u + 1] - uniq_digit_offsets[r][u];
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

      auto combine_simd = [&](int64_t d) -> cbatch {
        int const *info = map_ptr + d * map_stride;
        cbatch pow      = uq_simd_buf[0][info[0]];
        if (info[Rank] & 1) pow = xsimd::conj(pow);
        cbatch p1 = uq_simd_buf[1][info[1]];
        if (info[Rank] & 2) p1 = xsimd::conj(p1);
        pow *= p1;
        if constexpr (Rank > 2) {
          for (int r = 2; r < Rank; ++r) {
            cbatch pr = uq_simd_buf[r][info[r]];
            if (info[Rank] & (1 << r)) pr = xsimd::conj(pr);
            pow *= pr;
          }
        }
        return pow;
      };

      auto combine_scalar = [&](int64_t d) -> dcomplex {
        int const *info = map_ptr + d * map_stride;
        dcomplex pow    = uq_scalar_buf[0][info[0]];
        if (info[Rank] & 1) pow = std::conj(pow);
        for (int r = 1; r < Rank; ++r) {
          dcomplex pr = uq_scalar_buf[r][info[r]];
          if (info[Rank] & (1 << r)) pr = std::conj(pr);
          pow *= pr;
        }
        return pow;
      };

      auto accumulate_simd_range = [&](int j_begin, int j_end) {
        for (int j = j_begin; j < j_end; j += simd_size) {
          poet::static_for<Rank>([&](const auto r) {
            for (int u = 0; u < n_unique[r]; ++u) uq_simd_buf[r][u] = compute_unique_simd(r, u, j);
          });
          cbatch fj = cbatch::load_unaligned(state.fx_arr.data() + j);
          for (int64_t d = 0; d < state.n_targets; ++d) sums_buf[d] = xsimd::fma(fj, combine_simd(d), sums_buf[d]);
        }
      };

      auto reduce_sums = [&]() {
        for (int64_t d = 0; d < state.n_targets; ++d) fiw_ptr[d] += xsimd::reduce_add(sums_buf[d]);
      };

      if (use_blocking) {
        constexpr int source_block = 128;
        for (int jb = 0; jb < buf_counter_simd; jb += source_block) {
          std::fill(sums_buf.begin(), sums_buf.end(), cbatch(dcomplex{0, 0}));
          accumulate_simd_range(jb, std::min(jb + source_block, static_cast<int>(buf_counter_simd)));
          reduce_sums();
        }
      } else {
        std::fill(sums_buf.begin(), sums_buf.end(), cbatch(dcomplex{0, 0}));
        accumulate_simd_range(0, static_cast<int>(buf_counter_simd));
        reduce_sums();
      }

      // Scalar tail
      for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
        poet::static_for<Rank>([&](const auto r) {
          for (int u = 0; u < n_unique[r]; ++u) uq_scalar_buf[r][u] = compute_unique_scalar(r, u, j);
        });
        dcomplex fj = state.fx_arr[j];
        for (int64_t d = 0; d < state.n_targets; ++d) fiw_ptr[d] += fj * combine_scalar(d);
      }
    }
  };

} // namespace triqs::utility::nfft
