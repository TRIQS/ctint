// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <nda/nda.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include <bit>
#include <triqs/mesh/matsubara_freq.hpp>

#include "finufft.h"
#include <xsimd/xsimd.hpp>
#include <poet/poet.hpp>

namespace triqs::utility {

  inline void check_finufft(int err) {
    if (err > 0) NDA_RUNTIME_ERROR << "Error in FINUFFT: " << err << "\n";
  }

  using finufft_plan_ptr = std::unique_ptr<finufft_plan_s, decltype([](finufft_plan p) { if (p) finufft_destroy(p); })>;

  using nda::array_view;
  using dcomplex = std::complex<double>;

  enum class nfft_type_t { automatic, type1, type3, direct_type1, direct_type3, direct_bitwise, direct_prime };

  template <int Rank> struct nfft_buf_t {

    static_assert(Rank >= 1 and Rank <= 3, "nfft_buf_t only supports Rank 1, 2, and 3");

    nfft_buf_t() = default;

    /// Type 1: non-uniform tau -> uniform Matsubara grid
    nfft_buf_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, double tol_ = 1e-15)
       : fiw_arr(std::move(fiw_arr_)),
         niws(nda::stdutil::make_std_array<int64_t>(fiw_arr.shape())),
         buf_size(buf_size_),
         beta(beta_),
         x_arr(Rank, buf_size),
         fx_arr(buf_size),
         fk_arr(fiw_arr.shape()),
         tol(tol_) {

      for (int n : niws) {
        if (n % 2 != 0) NDA_RUNTIME_ERROR << " dimension with uneven frequency count not allowed in NFFT Buffer \n";
        common_factor *= (n / 2) % 2 ? -1 : 1;
      }

      finufft_default_opts(&opts);
      opts.nthreads = 1;
      auto Ns = std::vector(niws.rbegin(), niws.rend());
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(1, Rank, Ns.data(), 1, 1, tol, &raw_plan, &opts));
      plan.reset(raw_plan);
    }

    /// Non-uniform target constructor: automatic dispatch, FINUFFT type3, or direct DFT
    nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<std::array<mesh::matsubara_freq, Rank>> target_mf_, int buf_size_,
               nfft_type_t type = nfft_type_t::automatic, double tol_ = 1e-15)
       : nfft_type(type),
         fiw_vec(std::move(fiw_vec_)),
         buf_size(buf_size_),
         n_targets(static_cast<int64_t>(target_mf_.size())),
         x_arr(Rank, buf_size_),
         fx_arr(buf_size_),
         tol(tol_) {

      // Contiguous intermediate buffer for all non-uniform types.
      // Direct kernels write to fk_vec (contiguous), then copy to fiw_vec (possibly strided).
      fk_vec.resize(n_targets);

      if (type == nfft_type_t::type3) {
        init_type3(target_mf_);

      } else if (type == nfft_type_t::direct_type1) {
        init_direct_common(target_mf_);

        // Compute min/range per dimension for dense power table addressing
        for (int r = 0; r < Rank; ++r) {
          auto row       = target_n(r, nda::range::all);
          auto [mn, mx]  = std::ranges::minmax_element(row);
          n_min_arr[r]   = *mn;
          n_range_arr[r] = *mx - *mn + 1;
        }
        for (int r = 0; r < Rank; ++r) pow_tbl[r].resize(buf_size_, n_range_arr[r]);

        // Precompute doubled offset indices for AoS gather
        target_idx.resize(Rank * n_targets);
        for (int r = 0; r < Rank; ++r)
          for (int64_t d = 0; d < n_targets; ++d) target_idx[r * n_targets + d] = 2 * (target_n(r, d) - n_min_arr[r]);

      } else if (type == nfft_type_t::direct_type3) {
        init_direct_naf(target_mf_, buf_size_);

      } else if (type == nfft_type_t::direct_bitwise) {
        init_direct_common(target_mf_);
        unsigned long max_exponent = 0;
        for (int r = 0; r < Rank; ++r) {
          bitwise_pow2_bits[r].resize(n_targets);
          for (int64_t d = 0; d < n_targets; ++d) {
            unsigned long exponent = odd_exponent_abs(target_n(r, d));
            max_exponent           = std::max(max_exponent, exponent);
            std::vector<int> bits;
            for (int k = 0; exponent > 0; ++k, exponent >>= 1)
              if (exponent & 1ul) bits.push_back(k);
            bitwise_pow2_bits[r][d] = std::move(bits);
          }
        }
        num_pow2_levels_bitwise = std::max(1, static_cast<int>(std::bit_width(max_exponent)));
        for (int r = 0; r < Rank; ++r) bitwise_pow2_tbl[r].resize(num_pow2_levels_bitwise, buf_size_);

      } else if (type == nfft_type_t::direct_prime) {
        init_direct_common(target_mf_);
        std::vector<int> all_primes;
        for (int r = 0; r < Rank; ++r) {
          target_prime_sums[r].resize(n_targets);
          for (int64_t d = 0; d < n_targets; ++d) {
            target_prime_sums[r][d] = express_as_prime_sum(static_cast<long>(odd_exponent_abs(target_n(r, d))));
            for (int p : target_prime_sums[r][d]) all_primes.push_back(p);
          }
        }
        std::sort(all_primes.begin(), all_primes.end());
        all_primes.erase(std::unique(all_primes.begin(), all_primes.end()), all_primes.end());
        primes = std::move(all_primes);
        for (int r = 0; r < Rank; ++r)
          for (int64_t d = 0; d < n_targets; ++d)
            for (int &p : target_prime_sums[r][d]) p = static_cast<int>(std::find(primes.begin(), primes.end(), p) - primes.begin());
        for (int r = 0; r < Rank; ++r) prime_pow_tbl[r].resize(primes.size(), buf_size_);

      } else if (type == nfft_type_t::automatic) {
        init_type3(target_mf_);
        init_direct_naf(target_mf_, buf_size_);

      } else {
        NDA_RUNTIME_ERROR << "nfft_buf_t: unsupported nfft_type_t for non-uniform target constructor\n";
      }
    }

    /// Convenience constructor for Rank=1: accepts vector of matsubara_freq directly
    nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<mesh::matsubara_freq> const &target_mf_, int buf_size_,
               nfft_type_t type = nfft_type_t::automatic, double tol_ = 1e-15)
      requires(Rank == 1)
       : nfft_buf_t(
            std::move(fiw_vec_),
            [&] {
              std::vector<std::array<mesh::matsubara_freq, 1>> result;
              result.reserve(target_mf_.size());
              for (auto const &mf : target_mf_) result.push_back({mf});
              return result;
            }(),
            buf_size_, type, tol_) {}

    ~nfft_buf_t() {
      if (buf_counter != 0) std::cout << " WARNING: Points in NFFT Buffer lost \n";
    }

    nfft_buf_t(nfft_buf_t const &)            = delete;
    nfft_buf_t &operator=(nfft_buf_t const &) = delete;
    nfft_buf_t(nfft_buf_t &&)                 = default;
    nfft_buf_t &operator=(nfft_buf_t &&rhs) noexcept {
      // array_view::operator= does deep copy; we need rebind via destroy + placement new
      if (this != &rhs) {
        std::destroy_at(this);
        std::construct_at(this, std::move(rhs));
      }
      return *this;
    }

    void rebind(array_view<dcomplex, Rank> new_fiw_arr) {
      flush();
      TRIQS_ASSERT((new_fiw_arr.shape() == fiw_arr.shape() or fiw_arr.empty())
                   and " Nfft Buffer: Rebind to array of different shape not allowed ");
      fiw_arr.rebind(new_fiw_arr);
    }

    void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {
      if (x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      if (nfft_type == nfft_type_t::type1) {
        double tau_sum = 0.0;
        for (int r = 0; r < Rank; ++r) {
          x_arr(r, buf_counter) = 2 * M_PI * (tau_arr[r] / beta - 0.5);
          tau_sum += tau_arr[r];
        }
        fx_arr[buf_counter] = std::exp(dcomplex(0, M_PI * tau_sum / beta)) * ftau;
      } else {
        for (int r = 0; r < Rank; ++r) x_arr(r, buf_counter) = tau_arr[r];
        fx_arr[buf_counter] = ftau;
      }

      if (++buf_counter >= buf_size) {
        do_nfft();
        buf_counter = 0;
      }
    }

    void flush() {
      if (x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";
      if (buf_counter == 0) return;
      do_nfft();
      buf_counter = 0;
    }

    private:
    nfft_type_t nfft_type = nfft_type_t::type1;

    // Type 1 output (uniform grid)
    nda::array_view<dcomplex, Rank> fiw_arr;
    nda::array_view<dcomplex, 1> fiw_vec; // Type 3 / direct output (non-uniform targets)
    std::array<int64_t, Rank> niws{};
    finufft_plan_ptr plan;
    int buf_size = 0;
    double beta = 0;
    int buf_counter = 0;
    int common_factor = 1;
    int64_t n_targets = 0;
    finufft_opts opts{};
    nda::array<double, 2> x_arr;
    nda::vector<dcomplex> fx_arr;
    nda::array<dcomplex, Rank> fk_arr;
    nda::array<double, 2> s_arr;       // Type 3 target frequencies
    nda::vector<dcomplex> fk_vec;      // Type 3 output buffer
    nda::array<long, 2> target_n;      // Direct: integer Matsubara indices (Rank, n_targets)
    double tol = 1e-15;

    // direct_type1 (dense): per-dimension min index, range, power table, and offset indices
    std::array<long, Rank> n_min_arr{};
    std::array<long, Rank> n_range_arr{};
    mutable std::array<nda::array<dcomplex, 2>, Rank> pow_tbl; // (buf_size, n_range) per rank
    std::vector<long> target_idx;                               // doubled offset indices for AoS gather

    // Bitwise kernel: per-rank pow2 tables and bit lists
    int num_pow2_levels_bitwise = 0;
    std::array<nda::array<dcomplex, 2>, Rank> bitwise_pow2_tbl;
    std::array<std::vector<std::vector<int>>, Rank> bitwise_pow2_bits;

    // Prime-sum kernel: per-rank prime decomposition tables
    std::vector<int> primes;
    std::array<nda::array<dcomplex, 2>, Rank> prime_pow_tbl;
    std::array<std::vector<std::vector<int>>, Rank> target_prime_sums;

    // NAF kernel: per-rank pow2 table and flattened signed digit lists
    int naf_num_pow2_levels = 0;
    std::array<nda::array<dcomplex, 2>, Rank> naf_pow2_tbl;
    // Flattened per-rank digit data: all targets' digits stored contiguously.
    // Each digit is a row index; negative values encode conjugation: -(row+1).
    std::array<std::vector<int>, Rank> naf_digits_flat;
    // naf_digit_offsets[r][d] = start index in naf_digits_flat[r] for target d.
    // naf_digit_offsets[r][n_targets] = total number of digits.
    std::array<std::vector<int>, Rank> naf_digit_offsets;

    // Dispatch threshold for automatic mode: use direct_type3 (NAF) when
    // buf_counter * n_targets < threshold, otherwise fall back to FINUFFT type3.
    static constexpr int64_t dispatch_threshold = 50'000'000;

    using target_mf_vec = std::vector<std::array<mesh::matsubara_freq, Rank>>;

    void init_type3(target_mf_vec const &target_mf_) {
      s_arr.resize(Rank, n_targets);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < n_targets; ++d) s_arr(r, d) = std::imag(dcomplex(target_mf_[d][r]));
      finufft_default_opts(&opts);
      opts.nthreads         = 1;
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(3, Rank, nullptr, 1, 1, tol, &raw_plan, &opts));
      plan.reset(raw_plan);
    }

    void init_direct_common(target_mf_vec const &target_mf_) {
      beta = target_mf_[0][0].beta;
      target_n.resize(Rank, n_targets);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < n_targets; ++d) target_n(r, d) = target_mf_[d][r].n;
    }

    void init_direct_naf(target_mf_vec const &target_mf_, int buf_size_) {
      init_direct_common(target_mf_);
      unsigned long max_exponent = 0;
      for (int r = 0; r < Rank; ++r) {
        naf_digit_offsets[r].resize(n_targets + 1);
        naf_digit_offsets[r][0] = 0;
        for (int64_t d = 0; d < n_targets; ++d) {
          unsigned long exp = odd_exponent_abs(target_n(r, d));
          max_exponent      = std::max(max_exponent, exp);
          auto digits       = compute_naf(exp);
          naf_digits_flat[r].insert(naf_digits_flat[r].end(), digits.begin(), digits.end());
          naf_digit_offsets[r][d + 1] = static_cast<int>(naf_digits_flat[r].size());
        }
      }
      naf_num_pow2_levels = std::max(1, static_cast<int>(std::bit_width(max_exponent)) + 1);
      for (int r = 0; r < Rank; ++r) naf_pow2_tbl[r].resize(naf_num_pow2_levels, buf_size_);
    }

    // |2n+1| for fermionic Matsubara index n
    static constexpr unsigned long odd_exponent_abs(long n) {
      long odd = 2 * n + 1;
      return static_cast<unsigned long>(odd >= 0 ? odd : -odd);
    }

    static constexpr bool is_prime(long x) {
      if (x < 2) return false;
      if (x == 2) return true;
      if (x % 2 == 0) return false;
      for (long i = 3; i * i <= x; i += 2)
        if (x % i == 0) return false;
      return true;
    }

    using cbatch                        = xsimd::batch<dcomplex>;
    static constexpr std::size_t simd_size = cbatch::size;
    static constexpr int n_acc_bitwise     = 4; // ILP accumulators for Rank==1
    static constexpr int n_acc_prime       = 4; // ILP accumulators for Rank>1
    static constexpr int n_acc_naf         = 4; // ILP accumulators for NAF

    static std::vector<int> express_as_prime_sum(long n) {
      std::vector<int> out;
      while (n > 0) {
        if (n == 1) { out.push_back(1); break; }
        if (n <= 3) { out.push_back(static_cast<int>(n)); break; }
        if (n == 4) { out.push_back(2); out.push_back(2); break; }
        long p = n;
        while (p > 1 && !is_prime(p)) --p;
        out.push_back(static_cast<int>(p));
        n -= p;
      }
      return out;
    }

    // NAF (Non-Adjacent Form) decomposition of n into signed binary digits.
    // Returns encoded digits: k for +1 at bit k, -(k+1) for -1 at bit k.
    static std::vector<int> compute_naf(unsigned long n) {
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

    // SIMD+ILP target accumulation: processes n_acc targets simultaneously.
    // compute_simd_pow(d, j) returns SIMD batch of exp(i*omega_d*tau_j).
    // compute_scalar_pow(d, j) returns scalar version for the tail.
    template <int n_acc, typename SimdPowFunc, typename ScalarPowFunc>
    [[gnu::always_inline]] inline void accumulate_targets_ilp(int64_t n_targets_total, int64_t buf_counter_simd, dcomplex *fiw_ptr,
                                                               SimdPowFunc &&compute_simd_pow, ScalarPowFunc &&compute_scalar_pow) {
      auto accumulate_one = [&](int64_t d) {
        cbatch sum_vec(dcomplex{0, 0});
        for (int j = 0; j < buf_counter_simd; j += simd_size)
          sum_vec = xsimd::fma(cbatch::load_unaligned(fx_arr.data() + j), compute_simd_pow(d, j), sum_vec);
        dcomplex sum = xsimd::reduce_add(sum_vec);
        for (int j = buf_counter_simd; j < buf_counter; ++j)
          sum += fx_arr[j] * compute_scalar_pow(d, j);
        fiw_ptr[d] += sum;
      };

      int64_t const n_targets_main = (n_targets_total / n_acc) * n_acc;
      int64_t d                    = 0;

      for (; d < n_targets_main; d += n_acc) {
        std::array<cbatch, n_acc> sum_vecs;
        poet::static_for<n_acc>([&](const auto i) { sum_vecs[i] = cbatch(dcomplex{0, 0}); });

        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          cbatch fj = cbatch::load_unaligned(fx_arr.data() + j);
          poet::static_for<n_acc>([&](const auto i) {
            sum_vecs[i] = xsimd::fma(fj, compute_simd_pow(d + i, j), sum_vecs[i]);
          });
        }

        std::array<dcomplex, n_acc> sums;
        poet::static_for<n_acc>([&](const auto i) { sums[i] = xsimd::reduce_add(sum_vecs[i]); });

        for (int j = buf_counter_simd; j < buf_counter; ++j) {
          dcomplex fj = fx_arr[j];
          poet::static_for<n_acc>([&](const auto i) { sums[i] += fj * compute_scalar_pow(d + i, j); });
        }

        poet::static_for<n_acc>([&](const auto i) { fiw_ptr[d + i] += sums[i]; });
      }

      for (; d < n_targets_total; ++d) accumulate_one(d);
    }

    // Run a direct kernel into fk_vec (contiguous), then copy to fiw_vec (possibly strided).
    template <typename DirectKernel> void run_direct(DirectKernel &&kernel) {
      fk_vec = 0;
      kernel();
      fiw_vec += fk_vec;
    }

    void do_nfft() {
      if (nfft_type == nfft_type_t::type1)
        do_nfft_type1();
      else if (nfft_type == nfft_type_t::automatic) {
        if (static_cast<int64_t>(buf_counter) * n_targets < dispatch_threshold)
          run_direct([this] { do_direct_naf(); });
        else
          do_nfft_type3();
      } else if (nfft_type == nfft_type_t::type3)
        do_nfft_type3();
      else if (nfft_type == nfft_type_t::direct_type3)
        run_direct([this] { do_direct_naf(); });
      else if (nfft_type == nfft_type_t::direct_bitwise)
        run_direct([this] { do_direct_bitwise(); });
      else if (nfft_type == nfft_type_t::direct_prime)
        run_direct([this] { do_direct_prime(); });
      else
        run_direct([this] { do_direct_type1(); });
    }

    // FINUFFT expects coordinates in reverse rank order
    void set_pts(nda::array<double, 2> *tgt = nullptr) {
      auto _ = nda::range::all;
      auto n_tgt = tgt ? n_targets : int64_t{0};
      auto t = [&](int r) -> double * { return tgt ? (*tgt)(r, _).data() : nullptr; };
      if constexpr (Rank == 1)
        check_finufft(finufft_setpts(plan.get(), buf_counter, x_arr(0, _).data(), nullptr, nullptr, n_tgt, t(0), nullptr, nullptr));
      else if constexpr (Rank == 2)
        check_finufft(finufft_setpts(plan.get(), buf_counter, x_arr(1, _).data(), x_arr(0, _).data(), nullptr, n_tgt, t(1), t(0), nullptr));
      else
        check_finufft(finufft_setpts(plan.get(), buf_counter, x_arr(2, _).data(), x_arr(1, _).data(), x_arr(0, _).data(), n_tgt, t(2), t(1), t(0)));
    }

    void do_nfft_type1() {
      set_pts();
      check_finufft(finufft_execute(plan.get(), fx_arr.data(), fk_arr.data()));
      for (auto idx_tpl : fiw_arr.indices()) {
        auto idx_sum = std::apply([](auto... idx) { return (idx + ... + 0); }, idx_tpl);
        int factor   = common_factor * (idx_sum % 2 ? -1 : 1);
        std::apply(fiw_arr, idx_tpl) += std::apply(fk_arr, idx_tpl) * factor;
      }
    }

    void do_nfft_type3() {
      set_pts(&s_arr);
      check_finufft(finufft_execute(plan.get(), fx_arr.data(), fk_vec.data()));
      fiw_vec += fk_vec;
    }

    // Direct NUDFT via bitwise power-of-two decomposition.
    // z = exp(i*pi*tau/beta), z^|2n+1| computed via binary exponentiation.
    void do_direct_bitwise() {
      double const pi_over_beta      = M_PI / beta;
      int64_t const buf_counter_simd = buf_counter & -simd_size;
      int64_t const stride           = buf_size;
      dcomplex *fiw_ptr              = fk_vec.data();

      // Phase 1: Build per-rank pow2 tables via sincos + repeated squaring
      poet::static_for<Rank>([&](const auto r) {
        dcomplex *tbl = bitwise_pow2_tbl[r].data();

        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          using rbatch            = xsimd::batch<double>;
          auto [sin_vec, cos_vec] = xsimd::sincos(rbatch::load_unaligned(&x_arr(r, j)) * pi_over_beta);
          cbatch(cos_vec, sin_vec).store_unaligned(tbl + j);
        }
        for (int j = buf_counter_simd; j < buf_counter; ++j) {
          double const theta = pi_over_beta * x_arr(r, j);
          tbl[j]             = dcomplex{std::cos(theta), std::sin(theta)};
        }

        for (int k = 1; k < num_pow2_levels_bitwise; ++k) {
          dcomplex const *prev_row = tbl + (k - 1) * stride;
          dcomplex *cur_row        = tbl + k * stride;
          for (int j = 0; j < buf_counter_simd; j += simd_size) {
            cbatch prev = cbatch::load_unaligned(prev_row + j);
            (prev * prev).store_unaligned(cur_row + j);
          }
          for (int j = buf_counter_simd; j < buf_counter; ++j) {
            dcomplex prev = prev_row[j];
            cur_row[j]    = prev * prev;
          }
        }
      });

      // Phase 2: accumulate with per-rank bit products
      std::array<dcomplex const *, Rank> tbl_base;
      poet::static_for<Rank>([&](const auto r) { tbl_base[r] = bitwise_pow2_tbl[r].data(); });

      accumulate_targets_ilp<n_acc_bitwise>(
         n_targets, buf_counter_simd, fiw_ptr,
         [&](int64_t d, int j) -> cbatch {
           cbatch pow_prod;
           poet::static_for<Rank>([&](const auto r) {
             auto const *base = tbl_base[r];
             cbatch rank_pow(dcomplex{1.0, 0.0});
             for (int k : bitwise_pow2_bits[r][d]) rank_pow *= cbatch::load_unaligned(base + k * stride + j);
             rank_pow = target_n(r, d) < 0 ? xsimd::conj(rank_pow) : rank_pow;
             pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
           });
           return pow_prod;
         },
         [&](int64_t d, int j) -> dcomplex {
           dcomplex pow_prod{1.0, 0.0};
           poet::static_for<Rank>([&](const auto r) {
             auto const *base = tbl_base[r];
             dcomplex rank_pow{1.0, 0.0};
             for (int k : bitwise_pow2_bits[r][d]) rank_pow *= *(base + k * stride + j);
             rank_pow = target_n(r, d) < 0 ? std::conj(rank_pow) : rank_pow;
             pow_prod *= rank_pow;
           });
           return pow_prod;
         });
    }

    // Direct NUDFT via prime-sum decomposition.
    // |2n+1| = p1+p2+... so z^|2n+1| = z^p1 * z^p2 * ..., sharing z^p across targets.
    void do_direct_prime() {
      double const pi_over_beta      = M_PI / beta;
      int64_t const buf_counter_simd = buf_counter & -simd_size;
      dcomplex *fiw_ptr              = fk_vec.data();
      int const num_primes           = static_cast<int>(primes.size());

      // Build prime power table: prime_pow_tbl[r](p_idx, j) = z_r^prime
      poet::static_for<Rank>([&](const auto r) {
        std::vector<dcomplex> z_vals(buf_counter);
        for (int j = 0; j < buf_counter; ++j) {
          double const theta = pi_over_beta * x_arr(r, j);
          z_vals[j]          = dcomplex{std::cos(theta), std::sin(theta)};
        }

        for (int p_idx = 0; p_idx < num_primes; ++p_idx) {
          int prime = primes[p_idx];
          if (prime == 1) {
            for (int j = 0; j < buf_counter; ++j) prime_pow_tbl[r](p_idx, j) = z_vals[j];
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
            result.store_unaligned(&prime_pow_tbl[r](p_idx, j));
          }
          for (int j = buf_counter_simd; j < buf_counter; ++j) {
            dcomplex result{1.0, 0.0};
            dcomplex base = z_vals[j];
            for (int exp = prime; exp > 0; exp >>= 1) {
              if (exp & 1) result *= base;
              base *= base;
            }
            prime_pow_tbl[r](p_idx, j) = result;
          }
        }
      });

      accumulate_targets_ilp<n_acc_prime>(n_targets, buf_counter_simd, fiw_ptr,
        [&](int64_t d, int j) -> cbatch {
          cbatch pow_prod;
          poet::static_for<Rank>([&](const auto r) {
            cbatch rank_pow(dcomplex{1.0, 0.0});
            for (int pi : target_prime_sums[r][d]) rank_pow *= cbatch::load_unaligned(&prime_pow_tbl[r](pi, j));
            rank_pow = (target_n(r, d) < 0) ? xsimd::conj(rank_pow) : rank_pow;
            pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
          });
          return pow_prod;
        },
        [&](int64_t d, int j) -> dcomplex {
          dcomplex pow_prod{1.0, 0.0};
          poet::static_for<Rank>([&](const auto r) {
            dcomplex rank_pow{1.0, 0.0};
            for (int pi : target_prime_sums[r][d]) rank_pow *= prime_pow_tbl[r](pi, j);
            rank_pow  = (target_n(r, d) < 0) ? std::conj(rank_pow) : rank_pow;
            pow_prod *= rank_pow;
          });
          return pow_prod;
        });
    }

    // Rank-generic direct NUDFT via NAF (Non-Adjacent Form) decomposition.
    // Uses the same power-of-two table as bitwise, but signed digits {-1,0,+1}
    // reduce average non-zero count from L/2 to L/3. Conjugation for -1 digits is free.
    void do_direct_naf() {
      double const pi_over_beta      = M_PI / beta;
      int64_t const buf_counter_simd = buf_counter & -simd_size;
      int64_t const stride           = buf_size;
      dcomplex *fiw_ptr              = fk_vec.data();

      // Phase 1: Build per-rank pow2 tables via repeated squaring (raw pointer access)
      poet::static_for<Rank>([&](const auto r) {
        dcomplex *tbl = naf_pow2_tbl[r].data();

        // Level 0: z_r = exp(i*pi*tau_r/beta)
        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          using rbatch            = xsimd::batch<double>;
          auto [sin_vec, cos_vec] = xsimd::sincos(rbatch::load_unaligned(&x_arr(r, j)) * pi_over_beta);
          cbatch(cos_vec, sin_vec).store_unaligned(tbl + j);
        }
        for (int j = buf_counter_simd; j < buf_counter; ++j) {
          double const theta = pi_over_beta * x_arr(r, j);
          tbl[j]             = dcomplex{std::cos(theta), std::sin(theta)};
        }

        // Higher levels: squaring
        for (int k = 1; k < naf_num_pow2_levels; ++k) {
          dcomplex const *prev_row = tbl + (k - 1) * stride;
          dcomplex *cur_row        = tbl + k * stride;
          for (int j = 0; j < buf_counter_simd; j += simd_size) {
            cbatch prev = cbatch::load_unaligned(prev_row + j);
            (prev * prev).store_unaligned(cur_row + j);
          }
          for (int j = buf_counter_simd; j < buf_counter; ++j) {
            dcomplex prev = prev_row[j];
            cur_row[j]    = prev * prev;
          }
        }
      });

      // Phase 2: Source-blocked accumulation for cache locality.
      // The pow2 table is buf_size * n_levels * 16 bytes per rank, which can far exceed
      // L1/L2 cache. By processing sources in blocks, table data at each j position
      // stays in L1 as we iterate over all targets within the block.
      std::array<dcomplex const *, Rank> tbl_base;
      poet::static_for<Rank>([&](const auto r) { tbl_base[r] = naf_pow2_tbl[r].data(); });

      auto compute_simd_pow = [&](int64_t d, int j) -> cbatch {
        cbatch pow_prod;
        poet::static_for<Rank>([&](const auto r) {
          int const *digits = naf_digits_flat[r].data() + naf_digit_offsets[r][d];
          int n_digits      = naf_digit_offsets[r][d + 1] - naf_digit_offsets[r][d];
          auto const *base  = tbl_base[r];

          int d0          = digits[0];
          int row0        = d0 >= 0 ? d0 : -(d0 + 1);
          cbatch rank_pow = cbatch::load_unaligned(base + row0 * stride + j);
          if (d0 < 0) rank_pow = xsimd::conj(rank_pow);

          for (int i = 1; i < n_digits; ++i) {
            int di     = digits[i];
            int row    = di >= 0 ? di : -(di + 1);
            cbatch val = cbatch::load_unaligned(base + row * stride + j);
            rank_pow *= di >= 0 ? val : xsimd::conj(val);
          }

          rank_pow = (target_n(r, d) < 0) ? xsimd::conj(rank_pow) : rank_pow;
          pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
        });
        return pow_prod;
      };

      auto compute_scalar_pow = [&](int64_t d, int j) -> dcomplex {
        dcomplex pow_prod{1.0, 0.0};
        poet::static_for<Rank>([&](const auto r) {
          int const *digits = naf_digits_flat[r].data() + naf_digit_offsets[r][d];
          int n_digits      = naf_digit_offsets[r][d + 1] - naf_digit_offsets[r][d];
          auto const *base  = tbl_base[r];

          int d0            = digits[0];
          int row0          = d0 >= 0 ? d0 : -(d0 + 1);
          dcomplex rank_pow = *(base + row0 * stride + j);
          if (d0 < 0) rank_pow = std::conj(rank_pow);

          for (int i = 1; i < n_digits; ++i) {
            int di       = digits[i];
            int row      = di >= 0 ? di : -(di + 1);
            dcomplex val = *(base + row * stride + j);
            rank_pow *= di >= 0 ? val : std::conj(val);
          }

          rank_pow = (target_n(r, d) < 0) ? std::conj(rank_pow) : rank_pow;
          pow_prod *= rank_pow;
        });
        return pow_prod;
      };

      // Use source blocking when the table exceeds L2 cache, otherwise use unblocked ILP.
      constexpr int64_t l2_bytes           = 2 * 1024 * 1024;
      int64_t const table_bytes_per_source = Rank * naf_num_pow2_levels * static_cast<int64_t>(sizeof(dcomplex));
      bool const use_blocking              = buf_counter_simd * table_bytes_per_source > l2_bytes;

      if (use_blocking) {
        // Source-blocked: process sources in L1-sized blocks, iterating over all
        // targets per block so table data stays in L1 across target iterations.
        constexpr int source_block   = 128;
        constexpr int n_acc          = n_acc_naf;
        int64_t const n_targets_main = (n_targets / n_acc) * n_acc;

        for (int jb = 0; jb < buf_counter_simd; jb += source_block) {
          int const j_end = std::min(jb + source_block, static_cast<int>(buf_counter_simd));

          int64_t d = 0;
          for (; d < n_targets_main; d += n_acc) {
            std::array<cbatch, n_acc> local_sums;
            poet::static_for<n_acc>([&](const auto i) { local_sums[i] = cbatch(dcomplex{0, 0}); });

            for (int j = jb; j < j_end; j += simd_size) {
              cbatch fj = cbatch::load_unaligned(fx_arr.data() + j);
              poet::static_for<n_acc>([&](const auto i) { local_sums[i] = xsimd::fma(fj, compute_simd_pow(d + i, j), local_sums[i]); });
            }

            poet::static_for<n_acc>([&](const auto i) { fiw_ptr[d + i] += xsimd::reduce_add(local_sums[i]); });
          }
          for (; d < n_targets; ++d) {
            cbatch local_sum(dcomplex{0, 0});
            for (int j = jb; j < j_end; j += simd_size)
              local_sum = xsimd::fma(cbatch::load_unaligned(fx_arr.data() + j), compute_simd_pow(d, j), local_sum);
            fiw_ptr[d] += xsimd::reduce_add(local_sum);
          }
        }

        // Scalar tail: only needed for the blocking path (accumulate_targets_ilp handles it internally)
        for (int j = buf_counter_simd; j < buf_counter; ++j) {
          dcomplex fj = fx_arr[j];
          for (int64_t d = 0; d < n_targets; ++d) fiw_ptr[d] += fj * compute_scalar_pow(d, j);
        }
      } else {
        accumulate_targets_ilp<n_acc_naf>(n_targets, buf_counter_simd, fiw_ptr, compute_simd_pow, compute_scalar_pow);
      }
    }

    // Direct type1 (dense): geometric power table with SIMD gather
    void do_direct_type1() {
      using rbatch                       = xsimd::batch<double>;
      using ibatch                       = xsimd::batch<int64_t>;
      constexpr std::size_t simd_size_t1 = cbatch::size;

      double const pi_over_beta    = M_PI / beta;
      int64_t const n_targets_simd = n_targets - (n_targets % static_cast<int64_t>(simd_size_t1));
      dcomplex *fiw_ptr            = fk_vec.data();

      // Phase 1: Build power tables pow_tbl[r](j,i) = z^(2*(n_min+i)+1)
      for (int r = 0; r < Rank; ++r) {
        long const n_range      = n_range_arr[r];
        long const n_range_simd = n_range - (n_range % static_cast<long>(simd_size_t1));
        long const abs_n_min    = std::abs(n_min_arr[r]);
        bool const n_min_neg    = n_min_arr[r] < 0;

        for (int j = 0; j < buf_counter; ++j) {
          double const theta = pi_over_beta * x_arr(r, j);
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
      poet::static_for<0, Rank>([&](auto r) { idx_ptr[r] = target_idx.data() + r * n_targets; });

      for (int j = 0; j < buf_counter; ++j) {
        cbatch const fj(fx_arr[j]);

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
        for (int64_t d = n_targets_simd; d < n_targets; ++d) {
          dcomplex prod = 1.0;
          poet::static_for<0, Rank>([&](auto r) {
            prod *= reinterpret_cast<dcomplex const *>(pow_ptr[r])[idx_ptr[r][d] >> 1];
          });
          fiw_ptr[d] += fx_arr[j] * prod;
        }
      }
    }
  };

  template <int Rank>
  nfft_buf_t(nda::array_view<dcomplex, Rank>, int, double, double) -> nfft_buf_t<Rank>;

  template <std::size_t N>
  nfft_buf_t(nda::array_view<dcomplex, 1>, std::vector<std::array<mesh::matsubara_freq, N>>, int, nfft_type_t, double) -> nfft_buf_t<static_cast<int>(N)>;

  template <std::size_t N>
  nfft_buf_t(nda::array_view<dcomplex, 1>, std::vector<std::array<mesh::matsubara_freq, N>>, int) -> nfft_buf_t<static_cast<int>(N)>;

  nfft_buf_t(nda::array_view<dcomplex, 1>, std::vector<mesh::matsubara_freq> const &, int, nfft_type_t, double) -> nfft_buf_t<1>;

  nfft_buf_t(nda::array_view<dcomplex, 1>, std::vector<mesh::matsubara_freq> const &, int) -> nfft_buf_t<1>;

} // namespace triqs::utility
