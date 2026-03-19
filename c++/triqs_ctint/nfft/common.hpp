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
#include "sincos.hpp"
#include <xsimd/xsimd.hpp>
#include <poet/poet.hpp>

namespace triqs::utility::nfft {

  using nda::array_view;
  using dcomplex = std::complex<double>;

  template <int TolDigits = 12> inline dcomplex cis(double theta) {
    if constexpr (TolDigits >= 12) {
      double s, c;
      ::sincos(theta, &s, &c);
      return {c, s};
    } else {
      auto [s, c] = triqs::utility::math::sincos<TolDigits>(theta);
      return {c, s};
    }
  }

  inline void check_finufft(int err) {
    if (err > 0) NDA_RUNTIME_ERROR << "Error in FINUFFT: " << err << "\n";
  }

  using finufft_plan_ptr = std::unique_ptr<finufft_plan_s, decltype([](finufft_plan p) {
                             if (p) finufft_destroy(p);
                           })>;

  enum class type_t { automatic, type1, type1_gather, type3, direct_type1, direct_type3, direct_chain };

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

  using supported_tol_digits_t = std::integer_sequence<int, 6, 8, 10, 12>;

  constexpr int bucket_tol_digits(double tol) {
    if (tol <= 1e-12) return 12;
    if (tol <= 1e-10) return 10;
    if (tol <= 1e-8) return 8;
    return 6;
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

  // Transform raw tau coordinates in-place for FINUFFT type1 convention.
  // The hot part is embarrassingly parallel in the source index j, so run a SIMD
  // main loop over the contiguous x_arr rows and fx_arr, then finish with a scalar tail.
  template <int TolDigits = 12, int Rank> [[gnu::flatten]] void apply_type1_coord_transform(shared_state_t<Rank> &state, int n) {
    using dbatch = xsimd::batch<double>;
    using cbatch = xsimd::batch<dcomplex>;

    if (n <= 0) return;

    double const inv_beta = 1.0 / state.beta;
    constexpr int simd_width = static_cast<int>(dbatch::size);
    int const simd_n = n & -simd_width;

    std::array<double *, Rank> x_ptr{};
    for (int r = 0; r < Rank; ++r) x_ptr[r] = &state.x_arr(r, 0);
    dcomplex *fx_ptr = state.fx_arr.data();

    dbatch const inv_beta_vec(inv_beta);
    dbatch const half_vec(0.5);
    dbatch const two_pi_vec(2.0 * M_PI);
    dbatch const pi_over_beta_vec(M_PI * inv_beta);

    int j = 0;
    for (; j < simd_n; j += simd_width) {
      dbatch tau_sum(0.0);
      for (int r = 0; r < Rank; ++r) {
        dbatch tau = dbatch::load_unaligned(x_ptr[r] + j);
        tau_sum += tau;
        (two_pi_vec * (tau * inv_beta_vec - half_vec)).store_unaligned(x_ptr[r] + j);
      }

      auto [sin_theta, cos_theta] = triqs::utility::math::sincos<TolDigits>(pi_over_beta_vec * tau_sum);
      cbatch phase(cos_theta, sin_theta);
      (cbatch::load_unaligned(fx_ptr + j) * phase).store_unaligned(fx_ptr + j);
    }

    for (; j < n; ++j) {
      double tau_sum = 0.0;
      for (int r = 0; r < Rank; ++r) {
        double tau = state.x_arr(r, j);
        tau_sum += tau;
        state.x_arr(r, j) = 2 * M_PI * (tau * inv_beta - 0.5);
      }
      state.fx_arr[j] *= cis<TolDigits>(M_PI * tau_sum * inv_beta);
    }
  }

  // ---- SIMD type aliases and helpers ----

  using cbatch                           = xsimd::batch<dcomplex>;
  static constexpr std::size_t simd_size = cbatch::size;

  // ILP accumulator count: 2 for 16-reg ISAs (SSE/AVX2), 4 for 32-reg ISAs (AVX-512/NEON/SVE).
  static constexpr int n_acc = poet::vector_register_count() <= 16 ? 2 : 4;

  // SIMD+ILP target accumulation: processes n_acc targets simultaneously.
  // The inner source loop is shared across n_acc accumulators for ILP.
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
