// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../common.hpp"
#include <numeric>

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

      detect_sequential(state);

      if (is_sequential_) {
        // No power tables needed — fused path computes and accumulates in one pass
      } else {
        // Non-sequential path: allocate single-row power buffers and precompute gather indices
        for (int r = 0; r < Rank; ++r) pow_row[r].resize(n_range_arr[r]);

        if constexpr (Rank >= 2) {
          init_factored_groups(state);
        } else {
          target_idx.resize(state.n_targets);
          for (int64_t d = 0; d < state.n_targets; ++d) target_idx[d] = 2 * (state.target_n(0, d) - n_min_arr[0]);
        }
      }
    }

    void execute(shared_state_t<Rank> &state) {
      if (is_sequential_)
        execute_sequential(state);
      else
        execute_gather(state);
    }

    private:
    std::array<long, Rank> n_min_arr{};
    std::array<long, Rank> n_range_arr{};
    std::array<nda::vector<dcomplex>, Rank> pow_row;  // (n_range) per rank — buffer for gather path
    std::vector<long> target_idx;                      // doubled offset indices for AoS gather
    bool is_sequential_ = false;                        // true if targets form a dense Cartesian grid

    // Factored accumulation data (Rank >= 2): group targets by rank-0 index
    struct group_t {
      long pow0_offset;    // index into pow_row[0]
      int64_t start;       // start in padded arrays
      int64_t count;       // real targets in this group
      int64_t padded_count; // count rounded up to SIMD width
    };
    std::vector<group_t> groups_;
    std::vector<long> factored_idx_;   // permuted rank-1+ gather indices (doubled, AoS), padded
    std::vector<int64_t> scatter_map_; // scatter_map_[perm_d] = original target index (real entries only)
    int64_t n_padded_targets_ = 0;    // total padded target count
    nda::vector<dcomplex> grouped_fk_; // intermediate grouped output (padded size)

    // Build factored group data: sort targets by rank-0 index, group by unique n0 values.
    // Groups are padded to SIMD width to eliminate scalar tails in the hot loop.
    void init_factored_groups(shared_state_t<Rank> const &state) {
      constexpr int64_t S = static_cast<int64_t>(cbatch::size);

      // Sort target indices by rank-0 Matsubara index
      std::vector<int64_t> perm(state.n_targets);
      std::iota(perm.begin(), perm.end(), 0);
      std::sort(perm.begin(), perm.end(), [&](int64_t a, int64_t b) { return state.target_n(0, a) < state.target_n(0, b); });

      // Pass 1: identify groups and compute total padded size
      groups_.reserve(n_range_arr[0]);
      n_padded_targets_ = 0;
      int64_t g_start = 0;
      while (g_start < state.n_targets) {
        long n0       = state.target_n(0, perm[g_start]);
        int64_t g_end = g_start;
        while (g_end < state.n_targets && state.target_n(0, perm[g_end]) == n0) ++g_end;
        int64_t count        = g_end - g_start;
        int64_t padded_count = (count + S - 1) / S * S;
        groups_.push_back({n0 - n_min_arr[0], n_padded_targets_, count, padded_count});
        n_padded_targets_ += padded_count;
        g_start = g_end;
      }

      // Allocate arrays with padded sizes
      scatter_map_.resize(n_padded_targets_);
      factored_idx_.resize((Rank - 1) * n_padded_targets_);
      grouped_fk_.resize(n_padded_targets_);

      // Pass 2: fill gather indices and scatter map, with padding
      int64_t perm_pos = 0;
      for (auto const &g : groups_) {
        for (int64_t k = 0; k < g.count; ++k) {
          scatter_map_[g.start + k] = perm[perm_pos + k];
          for (int r = 1; r < Rank; ++r)
            factored_idx_[(r - 1) * n_padded_targets_ + g.start + k] = 2 * (state.target_n(r, perm[perm_pos + k]) - n_min_arr[r]);
        }
        // Pad with last real entry's index (result discarded by scatter)
        for (int64_t k = g.count; k < g.padded_count; ++k)
          for (int r = 1; r < Rank; ++r)
            factored_idx_[(r - 1) * n_padded_targets_ + g.start + k] = factored_idx_[(r - 1) * n_padded_targets_ + g.start + g.count - 1];
        perm_pos += g.count;
      }
    }

    // Detect if targets form a dense sequential grid (n_targets == product of n_range_arr)
    void detect_sequential(shared_state_t<Rank> const &state) {
      int64_t product = 1;
      for (int r = 0; r < Rank; ++r) product *= n_range_arr[r];
      if (product != state.n_targets) return;

      // Verify row-major ordering: target d maps to multi-index (k0, k1, ...) where
      // target_n(r, d) = n_min_arr[r] + k_r, with k_r cycling fastest for last rank
      for (int64_t d = 0; d < state.n_targets; ++d) {
        int64_t rem = d;
        for (int r = Rank - 1; r >= 0; --r) {
          long expected_n = n_min_arr[r] + static_cast<long>(rem % n_range_arr[r]);
          if (state.target_n(r, d) != expected_n) return;
          rem /= n_range_arr[r];
        }
      }
      is_sequential_ = true;
    }

    // Fast path: fused compute+accumulate for dense Cartesian grids.
    // Walks the geometric sequence exp(i*(2n+1)*pi*tau/beta) in-place,
    // accumulating directly into the output — no power table needed.
    void execute_sequential(shared_state_t<Rank> &state) {
      constexpr std::size_t S = cbatch::size;
      double const pi_over_beta = M_PI / state.beta;
      dcomplex * __restrict__ fiw_ptr = state.fk_vec.data();

      if constexpr (Rank == 1) {
        long const nr           = n_range_arr[0];
        long const nr_simd      = nr - (nr % static_cast<long>(S));
        double const base_freq  = static_cast<double>(2 * n_min_arr[0] + 1) * pi_over_beta;
        double const step_freq  = 2.0 * pi_over_beta;

        for (int j = 0; j < state.buf_counter; ++j) {
          double const tau        = state.x_arr(0, j);
          double const theta_base = base_freq * tau;
          double const theta_step = step_freq * tau;
          dcomplex const z_base = cis(theta_base);
          dcomplex const z_step = cis(theta_step);
          dcomplex const fj = state.fx_arr[j];

          // Build SIMD multiplier array: mult[k] = z_step^k
          alignas(cbatch::arch_type::alignment()) std::array<dcomplex, S> mult_arr;
          dcomplex z_pow = 1.0;
          for (std::size_t k = 0; k < S; ++k) { mult_arr[k] = z_pow; z_pow *= z_step; }
          cbatch const mult(cbatch::load_aligned(mult_arr.data()));
          cbatch const stride(z_pow); // z_step^S
          cbatch zp(fj * z_base);

          for (long d = 0; d < nr_simd; d += S) {
            xsimd::fma(zp, mult, cbatch::load_unaligned(fiw_ptr + d)).store_unaligned(fiw_ptr + d);
            zp *= stride;
          }
          // Scalar remainder
          dcomplex z_scalar = zp.get(0);
          for (long d = nr_simd; d < nr; ++d) { fiw_ptr[d] += z_scalar; z_scalar *= z_step; }
        }

      } else if constexpr (Rank == 2) {
        long const nr0          = n_range_arr[0];
        long const nr1          = n_range_arr[1];
        long const nr1_simd     = nr1 - (nr1 % static_cast<long>(S));
        double const base_freq0 = static_cast<double>(2 * n_min_arr[0] + 1) * pi_over_beta;
        double const step_freq0 = 2.0 * pi_over_beta;
        double const base_freq1 = static_cast<double>(2 * n_min_arr[1] + 1) * pi_over_beta;
        double const step_freq1 = 2.0 * pi_over_beta;

        for (int j = 0; j < state.buf_counter; ++j) {
          double const tau0 = state.x_arr(0, j);
          double const tau1 = state.x_arr(1, j);
          dcomplex const fj = state.fx_arr[j];

          double const th_base0 = base_freq0 * tau0;
          double const th_step0 = step_freq0 * tau0;
          dcomplex const z_base0 = cis(th_base0);
          dcomplex const z_step0 = cis(th_step0);

          double const th_base1 = base_freq1 * tau1;
          double const th_step1 = step_freq1 * tau1;
          dcomplex const z_base1 = cis(th_base1);
          dcomplex const z_step1 = cis(th_step1);

          // SIMD multiplier for inner dimension
          alignas(cbatch::arch_type::alignment()) std::array<dcomplex, S> mult_arr;
          dcomplex z_pow = 1.0;
          for (std::size_t k = 0; k < S; ++k) { mult_arr[k] = z_pow; z_pow *= z_step1; }
          cbatch const mult1(cbatch::load_aligned(mult_arr.data()));
          cbatch const stride1(z_pow); // z_step1^S

          dcomplex z0 = fj * z_base0;
          dcomplex * __restrict__ out = fiw_ptr;
          for (long k0 = 0; k0 < nr0; ++k0) {
            cbatch zp1(z0 * z_base1);
            for (long k1 = 0; k1 < nr1_simd; k1 += S) {
              xsimd::fma(zp1, mult1, cbatch::load_unaligned(out + k1)).store_unaligned(out + k1);
              zp1 *= stride1;
            }
            dcomplex z1_scalar = zp1.get(0);
            for (long k1 = nr1_simd; k1 < nr1; ++k1) { out[k1] += z1_scalar; z1_scalar *= z_step1; }
            out += nr1;
            z0 *= z_step0;
          }
        }

      } else {
        // Rank 3: fused triple nested loop
        long const nr0 = n_range_arr[0];
        long const nr1 = n_range_arr[1];
        long const nr2 = n_range_arr[2];
        long const nr2_simd = nr2 - (nr2 % static_cast<long>(S));
        double const base_freq0 = static_cast<double>(2 * n_min_arr[0] + 1) * pi_over_beta;
        double const step_freq0 = 2.0 * pi_over_beta;
        double const base_freq1 = static_cast<double>(2 * n_min_arr[1] + 1) * pi_over_beta;
        double const step_freq1 = 2.0 * pi_over_beta;
        double const base_freq2 = static_cast<double>(2 * n_min_arr[2] + 1) * pi_over_beta;
        double const step_freq2 = 2.0 * pi_over_beta;

        for (int j = 0; j < state.buf_counter; ++j) {
          double const tau0 = state.x_arr(0, j);
          double const tau1 = state.x_arr(1, j);
          double const tau2 = state.x_arr(2, j);
          dcomplex const fj = state.fx_arr[j];

          dcomplex const z_base0 = cis(base_freq0 * tau0);
          dcomplex const z_step0 = cis(step_freq0 * tau0);
          dcomplex const z_base1 = cis(base_freq1 * tau1);
          dcomplex const z_step1 = cis(step_freq1 * tau1);
          dcomplex const z_base2 = cis(base_freq2 * tau2);
          dcomplex const z_step2 = cis(step_freq2 * tau2);

          alignas(cbatch::arch_type::alignment()) std::array<dcomplex, S> mult_arr;
          dcomplex zp = 1.0;
          for (std::size_t k = 0; k < S; ++k) { mult_arr[k] = zp; zp *= z_step2; }
          cbatch const mult2(cbatch::load_aligned(mult_arr.data()));
          cbatch const stride2(zp);

          dcomplex z0 = fj * z_base0;
          dcomplex * __restrict__ out = fiw_ptr;
          for (long k0 = 0; k0 < nr0; ++k0) {
            dcomplex z01 = z0 * z_base1;
            for (long k1 = 0; k1 < nr1; ++k1) {
              cbatch zp2(z01 * z_base2);
              for (long k2 = 0; k2 < nr2_simd; k2 += S) {
                xsimd::fma(zp2, mult2, cbatch::load_unaligned(out + k2)).store_unaligned(out + k2);
                zp2 *= stride2;
              }
              dcomplex z2_scalar = zp2.get(0);
              for (long k2 = nr2_simd; k2 < nr2; ++k2) { out[k2] += z2_scalar; z2_scalar *= z_step2; }
              out += nr2;
              z01 *= z_step1;
            }
            z0 *= z_step0;
          }
        }
      }
    }

    // Compute z_step^n for |z_step|=1 using binary exponentiation + conjugate
    static dcomplex unit_pow(dcomplex z_step, long n) {
      if (n == 0) return {1.0, 0.0};
      bool neg            = (n < 0);
      unsigned long abs_n = static_cast<unsigned long>(neg ? -n : n);
      dcomplex result     = z_step;
      unsigned long p     = abs_n;
      // Find highest bit
      int bits = std::bit_width(p) - 1;
      for (int b = bits - 1; b >= 0; --b) {
        result *= result; // square
        if ((p >> b) & 1) result *= z_step;
      }
      return neg ? std::conj(result) : result;
    }

    // Fill destination with geometric sequence exp(i*(2n+1)*pi*tau/beta) for n in [n_min, n_min+n_range).
    void fill_pow_row(dcomplex * __restrict__ dest, int r, double tau, double pi_over_beta) {
      constexpr std::size_t S = cbatch::size;
      long const nr      = n_range_arr[r];
      long const nr_simd = nr - (nr % static_cast<long>(S));
      double const theta1 = pi_over_beta * tau;
      dcomplex const z1 = cis(theta1);
      dcomplex const z_step = z1 * z1;
      dcomplex const z_base = unit_pow(z_step, n_min_arr[r]) * z1;

      alignas(cbatch::arch_type::alignment()) std::array<dcomplex, S> mult_arr;
      dcomplex z_pow = 1.0;
      for (std::size_t k = 0; k < S; ++k) { mult_arr[k] = z_pow; z_pow *= z_step; }
      cbatch const mult(cbatch::load_aligned(mult_arr.data()));
      cbatch const stride(z_pow);
      cbatch zp_vec(z_base);

      for (long i = 0; i < nr_simd; i += S) {
        (zp_vec * mult).store_unaligned(dest + i);
        zp_vec *= stride;
      }
      dcomplex z_scalar = zp_vec.get(0);
      for (long i = nr_simd; i < nr; ++i) { dest[i] = z_scalar; z_scalar *= z_step; }
    }

    void build_pow_row(int r, double tau, double pi_over_beta) { fill_pow_row(pow_row[r].data(), r, tau, pi_over_beta); }

    // Rank=1: fused build+gather path with j-tiling to amortize output loads/stores.
    void execute_gather(shared_state_t<Rank> &state) requires(Rank == 1) {
      using rbatch              = xsimd::batch<double>;
      using ibatch              = xsimd::batch<int64_t>;
      constexpr std::size_t S   = cbatch::size;
      constexpr int B           = 4; // batch size for j-tiling
      double const pi_over_beta = M_PI / state.beta;

      int64_t const n_targets_simd = state.n_targets - (state.n_targets % static_cast<int64_t>(S));
      dcomplex * __restrict__ fiw_ptr = state.fk_vec.data();
      long const * __restrict__ idx0 = target_idx.data();

      std::array<nda::vector<dcomplex>, B> pow0_buf;
      for (int b = 0; b < B; ++b) pow0_buf[b].resize(n_range_arr[0]);

      int j = 0;
      for (; j + B <= state.buf_counter; j += B) {
        std::array<cbatch, B> fj_vec;
        for (int b = 0; b < B; ++b) {
          fill_pow_row(pow0_buf[b].data(), 0, state.x_arr(0, j + b), pi_over_beta);
          fj_vec[b] = cbatch(state.fx_arr[j + b]);
        }

        for (int64_t d = 0; d < n_targets_simd; d += S) {
          auto idx = ibatch::load_unaligned(idx0 + d);
          cbatch acc = cbatch::load_unaligned(fiw_ptr + d);
          for (int b = 0; b < B; ++b) {
            double const * __restrict__ p0 = reinterpret_cast<double const *>(pow0_buf[b].data());
            acc = xsimd::fma(fj_vec[b], cbatch(rbatch::gather(p0, idx), rbatch::gather(p0, idx + 1)), acc);
          }
          acc.store_unaligned(fiw_ptr + d);
        }
        for (int64_t d = n_targets_simd; d < state.n_targets; ++d)
          for (int b = 0; b < B; ++b)
            fiw_ptr[d] += state.fx_arr[j + b] * pow0_buf[b].data()[idx0[d] >> 1];
      }

      // Remainder
      for (; j < state.buf_counter; ++j) {
        fill_pow_row(pow0_buf[0].data(), 0, state.x_arr(0, j), pi_over_beta);
        cbatch const fj(state.fx_arr[j]);
        double const * __restrict__ pow0 = reinterpret_cast<double const *>(pow0_buf[0].data());

        for (int64_t d = 0; d < n_targets_simd; d += S) {
          auto idx = ibatch::load_unaligned(idx0 + d);
          cbatch pow_val(rbatch::gather(pow0, idx), rbatch::gather(pow0, idx + 1));
          xsimd::fma(fj, pow_val, cbatch::load_unaligned(fiw_ptr + d)).store_unaligned(fiw_ptr + d);
        }
        for (int64_t d = n_targets_simd; d < state.n_targets; ++d)
          fiw_ptr[d] += state.fx_arr[j] * pow0_buf[0].data()[idx0[d] >> 1];
      }
    }

    // Rank>=2: factored accumulation — group by rank-0 to eliminate one gather+multiply.
    // Groups are padded to SIMD width, eliminating scalar tails in the hot loop.
    // Buffer points are tiled in batches of B to amortize output loads/stores.
    void execute_gather(shared_state_t<Rank> &state) requires(Rank >= 2) {
      using rbatch              = xsimd::batch<double>;
      using ibatch              = xsimd::batch<int64_t>;
      constexpr std::size_t S   = cbatch::size;
      constexpr int B           = 4; // batch size for j-tiling
      double const pi_over_beta = M_PI / state.beta;

      grouped_fk_ = 0;
      dcomplex * __restrict__ gfk_ptr = grouped_fk_.data();

      if constexpr (Rank == 2) {
        // Allocate batch pow buffers (once per flush)
        std::array<nda::vector<dcomplex>, B> pow0_buf, pow1_buf;
        for (int b = 0; b < B; ++b) { pow0_buf[b].resize(n_range_arr[0]); pow1_buf[b].resize(n_range_arr[1]); }

        int j = 0;
        for (; j + B <= state.buf_counter; j += B) {
          // Build B sets of pow_rows
          for (int b = 0; b < B; ++b) {
            fill_pow_row(pow0_buf[b].data(), 0, state.x_arr(0, j + b), pi_over_beta);
            fill_pow_row(pow1_buf[b].data(), 1, state.x_arr(1, j + b), pi_over_beta);
          }

          for (auto const &g : groups_) {
            // Precompute B combined scalars: fj * pow0[group_offset]
            std::array<cbatch, B> combined_vec;
            for (int b = 0; b < B; ++b) combined_vec[b] = cbatch(state.fx_arr[j + b] * pow0_buf[b].data()[g.pow0_offset]);

            dcomplex * __restrict__ out    = gfk_ptr + g.start;
            long const * __restrict__ idx1 = factored_idx_.data() + g.start;

            for (int64_t k = 0; k < g.padded_count; k += S) {
              auto idx = ibatch::load_unaligned(idx1 + k);
              cbatch acc = cbatch::load_unaligned(out + k);
              // Accumulate B contributions with shared output load/store
              for (int b = 0; b < B; ++b) {
                double const * __restrict__ p1 = reinterpret_cast<double const *>(pow1_buf[b].data());
                acc = xsimd::fma(combined_vec[b], cbatch(rbatch::gather(p1, idx), rbatch::gather(p1, idx + 1)), acc);
              }
              acc.store_unaligned(out + k);
            }
          }
        }

        // Remainder: process one j at a time
        for (; j < state.buf_counter; ++j) {
          fill_pow_row(pow0_buf[0].data(), 0, state.x_arr(0, j), pi_over_beta);
          fill_pow_row(pow1_buf[0].data(), 1, state.x_arr(1, j), pi_over_beta);
          dcomplex const fj = state.fx_arr[j];

          for (auto const &g : groups_) {
            cbatch const combined_vec(fj * pow0_buf[0].data()[g.pow0_offset]);
            dcomplex * __restrict__ out    = gfk_ptr + g.start;
            double const * __restrict__ p1 = reinterpret_cast<double const *>(pow1_buf[0].data());
            long const * __restrict__ idx1 = factored_idx_.data() + g.start;

            for (int64_t k = 0; k < g.padded_count; k += S) {
              auto idx = ibatch::load_unaligned(idx1 + k);
              cbatch pv(rbatch::gather(p1, idx), rbatch::gather(p1, idx + 1));
              xsimd::fma(combined_vec, pv, cbatch::load_unaligned(out + k)).store_unaligned(out + k);
            }
          }
        }

      } else {
        // Rank > 2: original unbatched path
        for (int j = 0; j < state.buf_counter; ++j) {
          for (int r = 0; r < Rank; ++r) build_pow_row(r, state.x_arr(r, j), pi_over_beta);
          dcomplex const fj = state.fx_arr[j];

          for (auto const &g : groups_) {
            dcomplex const combined = fj * pow_row[0].data()[g.pow0_offset];
            cbatch const combined_vec(combined);
            dcomplex * __restrict__ out = gfk_ptr + g.start;

            for (int64_t k = 0; k < g.padded_count; k += S) {
              cbatch pow_prod;
              poet::static_for<1, Rank>([&](auto r) {
                double const *pr = reinterpret_cast<double const *>(pow_row[r].data());
                long const *idr  = factored_idx_.data() + (r - 1) * n_padded_targets_ + g.start;
                auto idx         = ibatch::load_unaligned(idr + k);
                cbatch pv(rbatch::gather(pr, idx), rbatch::gather(pr, idx + 1));
                pow_prod = (r == 1) ? pv : pow_prod * pv;
              });
              xsimd::fma(combined_vec, pow_prod, cbatch::load_unaligned(out + k)).store_unaligned(out + k);
            }
          }
        }
      }

      // Scatter only real (non-padded) entries to original target order
      dcomplex * __restrict__ fiw_ptr = state.fk_vec.data();
      for (auto const &g : groups_)
        for (int64_t k = 0; k < g.count; ++k) fiw_ptr[scatter_map_[g.start + k]] += gfk_ptr[g.start + k];
    }
  };

} // namespace triqs::utility::nfft
