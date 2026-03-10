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
      dcomplex *fiw_ptr              = state.fk_vec.data();
      int const num_primes           = static_cast<int>(primes.size());

      // Build prime power table: pow_tbl[r](p_idx, j) = z_r^prime
      poet::static_for<Rank>([&](const auto r) {
        std::vector<dcomplex> z_vals(state.buf_counter);
        for (int j = 0; j < state.buf_counter; ++j) {
          double const theta = pi_over_beta * state.x_arr(r, j);
          z_vals[j]          = dcomplex{std::cos(theta), std::sin(theta)};
        }

        for (int p_idx = 0; p_idx < num_primes; ++p_idx) {
          int prime = primes[p_idx];
          if (prime == 1) {
            for (int j = 0; j < state.buf_counter; ++j) pow_tbl[r](p_idx, j) = z_vals[j];
            continue;
          }
          // Binary exponentiation: z^prime
          for (int j = 0; j < buf_counter_simd; j += simd_size) {
            cbatch result(dcomplex{1.0, 0.0});
            cbatch base = cbatch::load_unaligned(&z_vals[j]);
            for (int exp = prime; exp > 0; exp >>= 1) {
              if (exp & 1) result *= base;
              base *= base;
            }
            result.store_unaligned(&pow_tbl[r](p_idx, j));
          }
          for (int j = buf_counter_simd; j < state.buf_counter; ++j) {
            dcomplex result{1.0, 0.0};
            dcomplex base = z_vals[j];
            for (int exp = prime; exp > 0; exp >>= 1) {
              if (exp & 1) result *= base;
              base *= base;
            }
            pow_tbl[r](p_idx, j) = result;
          }
        }
      });

      accumulate_targets_ilp<n_acc>(
         state.n_targets, buf_counter_simd, state.buf_counter, state.fx_arr.data(), fiw_ptr,
         [&](int64_t d, int j) -> cbatch {
           cbatch pow_prod;
           poet::static_for<Rank>([&](const auto r) {
             cbatch rank_pow(dcomplex{1.0, 0.0});
             for (int pi : prime_sums[r][d]) rank_pow *= cbatch::load_unaligned(&pow_tbl[r](pi, j));
             rank_pow = (state.target_n(r, d) < 0) ? xsimd::conj(rank_pow) : rank_pow;
             pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
           });
           return pow_prod;
         },
         [&](int64_t d, int j) -> dcomplex {
           dcomplex pow_prod{1.0, 0.0};
           poet::static_for<Rank>([&](const auto r) {
             dcomplex rank_pow{1.0, 0.0};
             for (int pi : prime_sums[r][d]) rank_pow *= pow_tbl[r](pi, j);
             rank_pow  = (state.target_n(r, d) < 0) ? std::conj(rank_pow) : rank_pow;
             pow_prod *= rank_pow;
           });
           return pow_prod;
         });
    }

    private:
    std::vector<int> primes;
    std::array<nda::array<dcomplex, 2>, Rank> pow_tbl;
    std::array<std::vector<std::vector<int>>, Rank> prime_sums;
  };

} // namespace triqs::utility::nfft
