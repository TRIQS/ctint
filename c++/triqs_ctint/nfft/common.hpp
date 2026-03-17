// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <nda/nda.hpp>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <memory>
#include <vector>
#include <triqs/mesh/matsubara_freq.hpp>

#include "finufft.h"
#include <xsimd/xsimd.hpp>
#include <poet/poet.hpp>

namespace triqs::utility::nfft {

  using nda::array_view;
  using dcomplex = std::complex<double>;

  // Combined sincos to avoid redundant trig computation (compiler doesn't merge without -ffast-math)
  inline dcomplex cis(double theta) {
    double s, c;
    ::sincos(theta, &s, &c);
    return {c, s};
  }

  inline void check_finufft(int err) {
    if (err > 0) NDA_RUNTIME_ERROR << "Error in FINUFFT: " << err << "\n";
  }

  using finufft_plan_ptr = std::unique_ptr<finufft_plan_s, decltype([](finufft_plan p) {
                             if (p) finufft_destroy(p);
                           })>;

  enum class type_t { automatic, type1, type1_gather, type3, direct_type1, direct_type3, direct_prime };

  using target_mf_t = mesh::matsubara_freq;

  // ---- Shared state between buffer_t and kernel structs ----

  template <int Rank> struct shared_state_t {
    int buf_size    = 0;
    double beta     = 0;
    int buf_counter = 0;
    nda::array<double, 2> x_arr;    // (Rank, buf_size)
    nda::vector<dcomplex> fx_arr;   // (buf_size)
    nda::vector<dcomplex> fk_vec;   // (n_targets) -- output for direct kernels
    int64_t n_targets = 0;
    nda::array<long, 2> target_n;   // (Rank, n_targets) -- Matsubara indices

    void init_direct_common(std::vector<std::array<target_mf_t, Rank>> const &target_mf) {
      beta = target_mf[0][0].beta;
      target_n.resize(Rank, n_targets);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < n_targets; ++d) target_n(r, d) = target_mf[d][r].n;
    }

    // Round up to SIMD boundary to eliminate scalar tail loops
    static constexpr int round_up_simd(int n) {
      constexpr int simd_sz = xsimd::batch<dcomplex>::size;
      return ((n + simd_sz - 1) / simd_sz) * simd_sz;
    }
  };

  // ---- Math utilities ----

  // |2n+1| for fermionic Matsubara index n
  constexpr unsigned long odd_exponent_abs(long n) {
    long odd = 2 * n + 1;
    return static_cast<unsigned long>(odd >= 0 ? odd : -odd);
  }

  constexpr bool is_prime(long x) {
    if (x < 2) return false;
    if (x == 2) return true;
    if (x % 2 == 0) return false;
    for (long i = 3; i * i <= x; i += 2)
      if (x % i == 0) return false;
    return true;
  }

  // Constexpr DP table for optimal (minimum-term) prime-sum decompositions.
  // best_summand[n] is the first summand to subtract; chase the chain to reconstruct.
  static constexpr int max_precomputed_exp = 2048;
  static constexpr auto prime_dp_table     = [] {
    std::array<uint16_t, max_precomputed_exp + 1> best{};
    std::array<uint16_t, max_precomputed_exp + 1> dp{};
    for (int i = 1; i <= max_precomputed_exp; ++i) dp[i] = 30000;

    // Sieve primes, collect into flat array
    std::array<bool, max_precomputed_exp + 1> sieve{};
    for (int i = 2; i <= max_precomputed_exp; ++i) sieve[i] = true;
    for (int i = 2; i * i <= max_precomputed_exp; ++i)
      if (sieve[i])
        for (int j = i * i; j <= max_precomputed_exp; j += i) sieve[j] = false;
    std::array<uint16_t, 320> plist{}; // pi(2048) = 309
    int np = 0;
    for (int i = 2; i <= max_precomputed_exp; ++i)
      if (sieve[i]) plist[np++] = static_cast<uint16_t>(i);

    // DP: dp[n] = minimum number of primes that sum to n
    for (int n = 1; n <= max_precomputed_exp; ++n) {
      if (dp[n - 1] + 1 < dp[n]) {
        dp[n]   = static_cast<uint16_t>(dp[n - 1] + 1);
        best[n] = 1;
      }
      for (int pi = 0; pi < np && plist[pi] <= n; ++pi) {
        int p = plist[pi];
        if (dp[n - p] + 1 < dp[n]) {
          dp[n]   = static_cast<uint16_t>(dp[n - p] + 1);
          best[n] = static_cast<uint16_t>(p);
        }
      }
    }
    return best;
  }();

  inline std::vector<int> express_as_prime_sum(long n) {
    std::vector<int> out;
    if (n <= max_precomputed_exp) {
      // Optimal decomposition from precomputed DP table
      while (n > 0) {
        int p = prime_dp_table[n];
        out.push_back(p);
        n -= p;
      }
    } else {
      // Greedy fallback for very large exponents
      while (n > 0) {
        if (n == 1) {
          out.push_back(1);
          break;
        }
        if (n <= 3) {
          out.push_back(static_cast<int>(n));
          break;
        }
        if (n == 4) {
          out.push_back(2);
          out.push_back(2);
          break;
        }
        long p = n;
        while (p > 1 && !is_prime(p)) --p;
        out.push_back(static_cast<int>(p));
        n -= p;
      }
    }
    return out;
  }

  // NAF (Non-Adjacent Form) decomposition of n into signed binary digits.
  // Returns encoded digits: k for +1 at bit k, -(k+1) for -1 at bit k.
  inline std::vector<int> compute_naf(unsigned long n) {
    std::vector<int> digits;
    long sn = static_cast<long>(n);
    for (int k = 0; sn > 0; ++k, sn >>= 1) {
      if (sn & 1) {
        int r = 2 - static_cast<int>(sn & 3); // +1 if sn%4==1, -1 if sn%4==3
        digits.push_back(r > 0 ? k : -(k + 1));
        sn -= r;
      }
    }
    return digits;
  }

  // Transform raw tau coordinates in-place for FINUFFT type1 convention
  template <int Rank> void apply_type1_coord_transform(shared_state_t<Rank> &state, int n) {
    double const inv_beta = 1.0 / state.beta;
    for (int j = 0; j < n; ++j) {
      double tau_sum = 0.0;
      for (int r = 0; r < Rank; ++r) {
        double tau = state.x_arr(r, j);
        tau_sum += tau;
        state.x_arr(r, j) = 2 * M_PI * (tau * inv_beta - 0.5);
      }
      state.fx_arr[j] *= cis(M_PI * tau_sum * inv_beta);
    }
  }

  // ---- SIMD type aliases and helpers ----

  using cbatch                           = xsimd::batch<dcomplex>;
  static constexpr std::size_t simd_size = cbatch::size;

  // Number of SIMD registers available on the target architecture.
  static constexpr int simd_num_regs = [] {
    using arch = typename cbatch::arch_type;
    // AVX-512 variants: 32 ZMM registers
    if constexpr (std::is_base_of_v<xsimd::avx512f, arch>) return 32;
#ifdef XSIMD_WITH_NEON
    // ARM NEON: 32 V registers
    else if constexpr (std::is_base_of_v<xsimd::neon, arch>) return 32;
#endif
#ifdef XSIMD_WITH_SVE
    // ARM SVE: 32 Z registers
    else if constexpr (std::is_base_of_v<xsimd::sve, arch>) return 32;
#endif
    // x86 SSE/AVX/AVX2: 16 XMM/YMM registers
    else return 16;
  }();

  // Optimal ILP accumulator count based on available SIMD registers.
  // Empirically determined from spill analysis (Clang 22, -O3):
  //   16 regs (SSE/AVX2): n_acc=2 is zero-spill; n_acc=3+ spills heavily
  //   32 regs (AVX-512/NEON/SVE): n_acc=4-5 is zero-spill
  static constexpr int n_acc = std::clamp(simd_num_regs <= 16 ? 2 : 4, 2, 8);

  // SIMD+ILP target accumulation: processes n_acc targets simultaneously.
  // compute_simd_pow(d, j) returns SIMD batch of exp(i*omega_d*tau_j).
  // compute_scalar_pow(d, j) returns scalar version for the tail.
  template <int n_acc, typename SimdPowFunc, typename ScalarPowFunc>
  [[gnu::always_inline]] inline void accumulate_targets_ilp(int64_t n_targets_total, int64_t buf_counter_simd, int buf_counter,
                                                            dcomplex *fx_data, dcomplex *fiw_ptr, SimdPowFunc &&compute_simd_pow,
                                                            ScalarPowFunc &&compute_scalar_pow) {
    auto accumulate_one = [&](int64_t d) {
      cbatch sum_vec(dcomplex{0, 0});
      for (int j = 0; j < buf_counter_simd; j += simd_size)
        sum_vec = xsimd::fma(cbatch::load_unaligned(fx_data + j), compute_simd_pow(d, j), sum_vec);
      dcomplex sum = xsimd::reduce_add(sum_vec);
      for (int j = buf_counter_simd; j < buf_counter; ++j) sum += fx_data[j] * compute_scalar_pow(d, j);
      fiw_ptr[d] += sum;
    };

    int64_t const n_targets_main = (n_targets_total / n_acc) * n_acc;
    int64_t d                    = 0;

    for (; d < n_targets_main; d += n_acc) {
      std::array<cbatch, n_acc> sum_vecs;
      poet::static_for<n_acc>([&](const auto i) { sum_vecs[i] = cbatch(dcomplex{0, 0}); });

      for (int j = 0; j < buf_counter_simd; j += simd_size) {
        cbatch fj = cbatch::load_unaligned(fx_data + j);
        poet::static_for<n_acc>([&](const auto i) { sum_vecs[i] = xsimd::fma(fj, compute_simd_pow(d + i, j), sum_vecs[i]); });
      }

      std::array<dcomplex, n_acc> sums;
      poet::static_for<n_acc>([&](const auto i) { sums[i] = xsimd::reduce_add(sum_vecs[i]); });

      for (int j = buf_counter_simd; j < buf_counter; ++j) {
        dcomplex fj = fx_data[j];
        poet::static_for<n_acc>([&](const auto i) { sums[i] += fj * compute_scalar_pow(d + i, j); });
      }

      poet::static_for<n_acc>([&](const auto i) { fiw_ptr[d + i] += sums[i]; });
    }

    for (; d < n_targets_total; ++d) accumulate_one(d);
  }

} // namespace triqs::utility::nfft
