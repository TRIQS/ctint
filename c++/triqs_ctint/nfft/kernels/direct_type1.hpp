// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../common.hpp"

namespace triqs::utility::nfft {

  template <int Rank> struct kernel_direct_type1_t {

    kernel_direct_type1_t() = default;

    kernel_direct_type1_t(shared_state_t<Rank> const &state, std::vector<std::array<target_mf_t, Rank>> const & /*target_mf*/) {
      // Compute min/range per dimension for dense power table addressing
      for (int r = 0; r < Rank; ++r) {
        auto row      = state.target_n(r, nda::range::all);
        auto [mn, mx] = std::ranges::minmax_element(row);
        n_min_arr[r]  = *mn;
        n_range_arr[r] = *mx - *mn + 1;
      }
      for (int r = 0; r < Rank; ++r) pow_tbl[r].resize(state.buf_size, n_range_arr[r]);

      // Precompute doubled offset indices for AoS gather
      target_idx.resize(Rank * state.n_targets);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < state.n_targets; ++d) target_idx[r * state.n_targets + d] = 2 * (state.target_n(r, d) - n_min_arr[r]);
    }

    void execute(shared_state_t<Rank> &state) {
      using rbatch                       = xsimd::batch<double>;
      using ibatch                       = xsimd::batch<int64_t>;
      constexpr std::size_t simd_size_t1 = cbatch::size;

      double const pi_over_beta    = M_PI / state.beta;
      int64_t const n_targets_simd = state.n_targets - (state.n_targets % static_cast<int64_t>(simd_size_t1));
      dcomplex *fiw_ptr            = state.fk_vec.data();

      // Phase 1: Build power tables pow_tbl[r](j,i) = z^(2*(n_min+i)+1)
      for (int r = 0; r < Rank; ++r) {
        long const n_range      = n_range_arr[r];
        long const n_range_simd = n_range - (n_range % static_cast<long>(simd_size_t1));
        long const abs_n_min    = std::abs(n_min_arr[r]);
        bool const n_min_neg    = n_min_arr[r] < 0;

        for (int j = 0; j < state.buf_counter; ++j) {
          double const theta = pi_over_beta * state.x_arr(r, j);
          dcomplex const z{std::cos(theta), std::sin(theta)};
          dcomplex const z2 = z * z;

          // Compute z^(2*n_min+1) via binary exponentiation
          dcomplex zp = z, base = n_min_neg ? std::conj(z2) : z2;
          for (long exp = abs_n_min; exp > 0; exp >>= 1) {
            if (exp & 1) zp *= base;
            base *= base;
          }

          // Fill geometric sequence: pow_row[i] = zp * z2^i using SIMD
          dcomplex *pow_row = &pow_tbl[r](j, 0);
          alignas(cbatch::arch_type::alignment()) std::array<dcomplex, simd_size_t1> mult_arr;
          dcomplex z2_pow = 1.0;
          for (std::size_t k = 0; k < simd_size_t1; ++k) {
            mult_arr[k] = z2_pow;
            z2_pow *= z2;
          }
          cbatch const mult(cbatch::load_aligned(mult_arr.data()));
          cbatch const stride(z2_pow);
          cbatch zp_vec(zp);

          for (long i = 0; i < n_range_simd; i += simd_size_t1) {
            (zp_vec * mult).store_unaligned(pow_row + i);
            zp_vec *= stride;
          }
          for (long i = n_range_simd; i < n_range; ++i) {
            pow_row[i] = zp_vec.get(0);
            zp_vec *= cbatch(z2);
          }
        }
      }

      // Phase 2: Accumulate f(tau) * product_r(pow_tbl[r][j][idx[r][d]]) into output
      std::array<long const *, Rank> idx_ptr;
      poet::static_for<0, Rank>([&](auto r) { idx_ptr[r] = target_idx.data() + r * state.n_targets; });

      for (int j = 0; j < state.buf_counter; ++j) {
        cbatch const fj(state.fx_arr[j]);

        std::array<double const *, Rank> pow_ptr;
        poet::static_for<0, Rank>([&](auto r) { pow_ptr[r] = reinterpret_cast<double const *>(&pow_tbl[r](j, 0)); });

        // SIMD loop with gather
        for (int64_t d = 0; d < n_targets_simd; d += simd_size_t1) {
          cbatch pow_prod;
          poet::static_for<0, Rank>([&](auto r) {
            auto idx = ibatch::load_unaligned(idx_ptr[r] + d);
            cbatch pow_val(rbatch::gather(pow_ptr[r], idx), rbatch::gather(pow_ptr[r], idx + 1));
            pow_prod = (r == 0) ? pow_val : pow_prod * pow_val;
          });
          xsimd::fma(fj, pow_prod, cbatch::load_unaligned(fiw_ptr + d)).store_unaligned(fiw_ptr + d);
        }

        // Scalar remainder
        for (int64_t d = n_targets_simd; d < state.n_targets; ++d) {
          dcomplex prod = 1.0;
          poet::static_for<0, Rank>(
             [&](auto r) { prod *= reinterpret_cast<dcomplex const *>(pow_ptr[r])[idx_ptr[r][d] >> 1]; });
          fiw_ptr[d] += state.fx_arr[j] * prod;
        }
      }
    }

    private:
    std::array<long, Rank> n_min_arr{};
    std::array<long, Rank> n_range_arr{};
    std::array<nda::array<dcomplex, 2>, Rank> pow_tbl; // (buf_size, n_range) per rank
    std::vector<long> target_idx;                       // doubled offset indices for AoS gather
  };

} // namespace triqs::utility::nfft
