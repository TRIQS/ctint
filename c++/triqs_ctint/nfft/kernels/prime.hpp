// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../common.hpp"

namespace triqs::utility::nfft {

  template <int Rank> struct kernel_prime_t {

    kernel_prime_t() = default;

    kernel_prime_t(shared_state_t<Rank> const &state) {
      std::vector<int> all_primes;
      for (int r = 0; r < Rank; ++r) {
        prime_sums[r].resize(state.n_targets);
        for (int64_t d = 0; d < state.n_targets; ++d) {
          prime_sums[r][d] = express_as_prime_sum(static_cast<long>(odd_exponent_abs(state.target_n(r, d))));
          for (int p : prime_sums[r][d]) all_primes.push_back(p);
        }
      }
      std::sort(all_primes.begin(), all_primes.end());
      all_primes.erase(std::unique(all_primes.begin(), all_primes.end()), all_primes.end());
      primes = std::move(all_primes);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < state.n_targets; ++d)
          for (int &p : prime_sums[r][d])
            p = static_cast<int>(std::find(primes.begin(), primes.end(), p) - primes.begin());
      for (int r = 0; r < Rank; ++r) pow_tbl[r].resize(primes.size(), state.buf_size);
    }

    void execute(shared_state_t<Rank> &state) {
      double const pi_over_beta      = M_PI / state.beta;
      int64_t const buf_counter_simd = state.buf_counter & -simd_size;
      int const buf_counter_padded   = std::min(shared_state_t<Rank>::round_up_simd(state.buf_counter), state.buf_size);
      dcomplex *fiw_ptr              = state.fk_vec.data();
      int const num_primes           = static_cast<int>(primes.size());
      int64_t const stride           = state.buf_size;

      // Build prime power table: pow_tbl[r](p_idx, j) = z_r^prime
      // Chain from previous prime: z^primes[i] = z^primes[i-1] * z^gap, where gap
      // is typically 2 (twin primes), so most entries cost just 2-3 multiplies.
      poet::static_for<Rank>([&](const auto r) {
        dcomplex *tbl = pow_tbl[r].data();

        // Row 0: z = exp(i*pi*tau/beta) via SIMD sincos
        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          using rbatch            = xsimd::batch<double>;
          auto [sin_vec, cos_vec] = xsimd::sincos(rbatch::load_unaligned(&state.x_arr(r, j)) * pi_over_beta);
          cbatch(cos_vec, sin_vec).store_unaligned(tbl + j);
        }
        // Scalar tail for row 0
        for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
          double const theta = pi_over_beta * state.x_arr(r, j);
          tbl[j]             = dcomplex{std::cos(theta), std::sin(theta)};
        }

        // Build z^prime for each prime by chaining from the previous prime's row.
        // When primes[0]==1, row 0 already holds z = z^1, so we start from p_idx=1.
        // z^primes[i] = z^primes[i-1] * z^(primes[i] - primes[i-1])
        for (int p_idx = (primes[0] == 1) ? 1 : 0; p_idx < num_primes; ++p_idx) {
          int prime = primes[p_idx];
          if (p_idx == 0 || primes[p_idx - 1] == 1) {
            // No useful predecessor: binary exp from z (row 0)
            for (int j = 0; j < buf_counter_simd; j += simd_size) {
              cbatch result(dcomplex{1.0, 0.0});
              cbatch base = cbatch::load_unaligned(tbl + j);
              for (int exp = prime; exp > 0; exp >>= 1) {
                if (exp & 1) result *= base;
                base *= base;
              }
              result.store_unaligned(tbl + p_idx * stride + j);
            }
            for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
              dcomplex result{1.0, 0.0};
              dcomplex base = tbl[j];
              for (int exp = prime; exp > 0; exp >>= 1) {
                if (exp & 1) result *= base;
                base *= base;
              }
              tbl[p_idx * stride + j] = result;
            }
          } else {
            int gap                    = prime - primes[p_idx - 1];
            dcomplex const *prev_row = tbl + (p_idx - 1) * stride;
            if (gap == 1) {
              // z^p = z^(p-1) * z (one multiply)
              for (int j = 0; j < buf_counter_simd; j += simd_size)
                (cbatch::load_unaligned(prev_row + j) * cbatch::load_unaligned(tbl + j)).store_unaligned(tbl + p_idx * stride + j);
              for (int j = buf_counter_simd; j < state.buf_counter; ++j) tbl[p_idx * stride + j] = prev_row[j] * tbl[j];
            } else if (gap == 2) {
              // z^p = z^(p-2) * z^2 (twin primes - most common case)
              for (int j = 0; j < buf_counter_simd; j += simd_size) {
                cbatch z = cbatch::load_unaligned(tbl + j);
                (cbatch::load_unaligned(prev_row + j) * z * z).store_unaligned(tbl + p_idx * stride + j);
              }
              for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
                dcomplex z                = tbl[j];
                tbl[p_idx * stride + j] = prev_row[j] * z * z;
              }
            } else {
              // General gap: binary exp of gap, multiply with prev
              for (int j = 0; j < buf_counter_simd; j += simd_size) {
                cbatch result(dcomplex{1.0, 0.0});
                cbatch base = cbatch::load_unaligned(tbl + j);
                for (int exp = gap; exp > 0; exp >>= 1) {
                  if (exp & 1) result *= base;
                  base *= base;
                }
                (cbatch::load_unaligned(prev_row + j) * result).store_unaligned(tbl + p_idx * stride + j);
              }
              for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
                dcomplex result{1.0, 0.0};
                dcomplex base = tbl[j];
                for (int exp = gap; exp > 0; exp >>= 1) {
                  if (exp & 1) result *= base;
                  base *= base;
                }
                tbl[p_idx * stride + j] = prev_row[j] * result;
              }
            }
          }
        }
      });

      // Accumulation with source blocking for L2 cache efficiency
      std::array<dcomplex const *, Rank> tbl_base;
      poet::static_for<Rank>([&](const auto r) { tbl_base[r] = pow_tbl[r].data(); });

      auto compute_simd_pow = [&](int64_t d, int j) -> cbatch {
        cbatch pow_prod;
        poet::static_for<Rank>([&](const auto r) {
          auto const *base    = tbl_base[r];
          auto const &digits  = prime_sums[r][d];
          int const n_digits = static_cast<int>(digits.size());

          // Skip identity initialization: load first digit directly
          cbatch rank_pow = cbatch::load_unaligned(base + digits[0] * stride + j);
          for (int i = 1; i < n_digits; ++i) rank_pow *= cbatch::load_unaligned(base + digits[i] * stride + j);

          rank_pow = (state.target_n(r, d) < 0) ? xsimd::conj(rank_pow) : rank_pow;
          pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
        });
        return pow_prod;
      };

      auto compute_scalar_pow = [&](int64_t d, int j) -> dcomplex {
        dcomplex pow_prod{1.0, 0.0};
        poet::static_for<Rank>([&](const auto r) {
          auto const &digits = prime_sums[r][d];
          int const n_digits = static_cast<int>(digits.size());

          dcomplex rank_pow = pow_tbl[r](digits[0], j);
          for (int i = 1; i < n_digits; ++i) rank_pow *= pow_tbl[r](digits[i], j);

          rank_pow  = (state.target_n(r, d) < 0) ? std::conj(rank_pow) : rank_pow;
          pow_prod *= rank_pow;
        });
        return pow_prod;
      };

      // L2-aware source blocking
      constexpr int64_t l2_bytes           = 2 * 1024 * 1024;
      int64_t const table_bytes_per_source = Rank * num_primes * static_cast<int64_t>(sizeof(dcomplex));
      bool const use_blocking              = buf_counter_padded * table_bytes_per_source > l2_bytes;

      if (use_blocking) {
        constexpr int source_block = 128;
        for (int jb = 0; jb < buf_counter_padded; jb += source_block) {
          int const j_end = std::min(jb + source_block, buf_counter_padded);
          // Accumulate for this block
          for (int64_t d = 0; d < state.n_targets; ++d) {
            cbatch sum_vec(dcomplex{0, 0});
            for (int j = jb; j < j_end && j < buf_counter_simd; j += simd_size)
              sum_vec = xsimd::fma(cbatch::load_unaligned(state.fx_arr.data() + j), compute_simd_pow(d, j), sum_vec);
            dcomplex sum = xsimd::reduce_add(sum_vec);
            for (int j = std::max(jb, static_cast<int>(buf_counter_simd)); j < std::min(j_end, state.buf_counter); ++j)
              sum += state.fx_arr[j] * compute_scalar_pow(d, j);
            fiw_ptr[d] += sum;
          }
        }
      } else {
        accumulate_targets_ilp<n_acc>(state.n_targets, buf_counter_simd, state.buf_counter, state.fx_arr.data(), fiw_ptr,
                                      compute_simd_pow, compute_scalar_pow);
      }
    }

    private:
    std::vector<int> primes;
    std::array<nda::array<dcomplex, 2>, Rank> pow_tbl;
    std::array<std::vector<std::vector<int>>, Rank> prime_sums;
  };

} // namespace triqs::utility::nfft
