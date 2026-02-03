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
#include <span>
#include <vector>

#include <triqs/mesh/matsubara_freq.hpp>

#include "finufft.h"
#include <xsimd/xsimd.hpp>
#include <poet/poet.hpp>

namespace triqs::utility {

  inline void check_finufft(int err) {
    if (err > 0) NDA_RUNTIME_ERROR << "Error in FINUFFT: " << err << "\n";
  }

  // RAII wrapper for finufft_plan (which is finufft_plan_s*)
  using finufft_plan_ptr = std::unique_ptr<finufft_plan_s, decltype([](finufft_plan p) { if (p) finufft_destroy(p); })>;

  using nda::array_view;
  using dcomplex = std::complex<double>;

  enum class nfft_type_t { type1, type3, direct };

  // Helper to convert vector<T> to vector<array<T, 1>> for Rank=1 convenience constructor
  inline std::vector<std::array<mesh::matsubara_freq, 1>> to_array_vector(std::vector<mesh::matsubara_freq> const &v) {
    std::vector<std::array<mesh::matsubara_freq, 1>> result;
    result.reserve(v.size());
    for (auto const &mf : v) result.push_back({mf});
    return result;
  }

  template <int Rank> struct nfft_buf_t {

    static_assert(Rank >= 1 and Rank <= 3, "nfft_buf_t only supports Rank 1, 2, and 3");

    /// Default constructor, creates unusable buffer!
    nfft_buf_t() = default;

    /**
     * Type 1 constructor: non-uniform tau -> uniform Matsubara grid
     *
     * Transforms scattered imaginary-time points to a regular Matsubara frequency grid.
     * Use this when accumulating into a full frequency mesh (e.g., G(iω_n) or G(iω_n, iω_m)).
     *
     * @param fiw_arr_  Output array of shape (N_1, ..., N_Rank) representing the uniform
     *                  Matsubara grid. Each dimension N_r must be even (fermionic frequencies).
     *                  Results are accumulated into this array.
     * @param buf_size_ Number of (tau, f(tau)) points to buffer before executing transform.
     * @param beta_     Inverse temperature. All tau values must be in [0, beta).
     * @param tol_      FINUFFT tolerance (default 1e-15).
     */
    nfft_buf_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, double tol_ = 1e-15)
       : fiw_arr(std::move(fiw_arr_)),
         niws(nda::stdutil::make_std_array<int64_t>(fiw_arr.shape())),
         buf_size(buf_size_),
         beta(beta_),
         x_arr(Rank, buf_size),
         fx_arr(buf_size),
         fk_arr(fiw_arr.shape()),
         tol(tol_) {

      // Capture frequency extents from fiw_arr and check that they are even ( i.e. fermionic matsubaras )
      for (int n : niws) {
        if (n % 2 != 0) NDA_RUNTIME_ERROR << " dimension with uneven frequency count not allowed in NFFT Buffer \n";
        common_factor *= (n / 2) % 2 ? -1 : 1; // Additional Minus sign for uneven Matsubara offset
      }

      // Init nfft_plan
      finufft_default_opts(&opts); // set default opts (must start with this)
      opts.nthreads = 1;           // enforce single-thread
      auto Ns = std::vector(niws.rbegin(), niws.rend()); // Reverse order for FINUFFT
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(/*type =*/1, Rank, Ns.data(), /*iflag=*/1, /*n_transf =*/1, tol, &raw_plan, &opts));
      plan.reset(raw_plan);
    }

    /**
     * Non-uniform target constructor: type3 (FINUFFT) or direct DFT
     *
     * Transforms scattered imaginary-time points to arbitrary (non-uniform) Matsubara frequencies.
     * Use this when you only need specific frequency points (e.g., DLR nodes, sparse sampling).
     *
     * @param fiw_vec_    Output vector of length n_targets. Results are accumulated here.
     * @param target_mf_  Vector of target frequency points. Each entry is an array of Rank
     *                    matsubara_freq objects specifying one multi-dimensional frequency point.
     *                    For Rank=2: target_mf_[d] = {iω_n, iν_m} gives the d-th target point.
     * @param buf_size_   Number of (tau, f(tau)) points to buffer before executing transform.
     * @param type        Transform algorithm:
     *                    - nfft_type_t::type3: Use FINUFFT type 3 (good for many targets)
     *                    - nfft_type_t::direct: Explicit DFT with SIMD (good for few targets)
     * @param tol_        FINUFFT tolerance for type3 (default 1e-15). Ignored for direct.
     */
    nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<std::array<mesh::matsubara_freq, Rank>> target_mf_, int buf_size_,
               nfft_type_t type, double tol_ = 1e-15)
       : nfft_type(type),
         fiw_vec(std::move(fiw_vec_)),
         buf_size(buf_size_),
         n_targets(static_cast<int64_t>(target_mf_.size())),
         x_arr(Rank, buf_size_),
         fx_arr(buf_size_),
         tol(tol_) {

      if (type == nfft_type_t::type3) {
        // Extract frequencies from matsubara_freq for FINUFFT type 3
        s_arr.resize(Rank, n_targets);
        for (int r = 0; r < Rank; ++r)
          for (int64_t d = 0; d < n_targets; ++d) s_arr(r, d) = std::imag(dcomplex(target_mf_[d][r]));
        fk_vec.resize(n_targets);
        finufft_default_opts(&opts);
        opts.nthreads = 1;
        finufft_plan raw_plan = nullptr;
        check_finufft(finufft_makeplan(3, Rank, nullptr, /*iflag=*/1, /*n_transf=*/1, tol, &raw_plan, &opts));
        plan.reset(raw_plan);

      } else if (type == nfft_type_t::direct) {
        // Extract integer indices and beta from matsubara_freq
        beta = target_mf_[0][0].beta;
        target_n.resize(Rank, n_targets);
        for (int r = 0; r < Rank; ++r)
          for (int64_t d = 0; d < n_targets; ++d) target_n(r, d) = target_mf_[d][r].n;
        // Compute min/range per dimension for power table addressing
        for (int r = 0; r < Rank; ++r) {
          auto row       = target_n(r, nda::range::all);
          auto [mn, mx]  = std::ranges::minmax_element(row);
          n_min_arr[r]   = *mn;
          n_range_arr[r] = *mx - *mn + 1;
        }
        // Preallocate power tables: shape (buf_size, n_range_arr[r]) per rank
        for (int r = 0; r < Rank; ++r) pow_tbl[r].resize(buf_size_, n_range_arr[r]);

        // Precompute offset indices (doubled): target_idx[r * n_targets + d] = 2 * (target_n(r, d) - n_min_arr[r])
        target_idx.resize(Rank * n_targets);
        for (int r = 0; r < Rank; ++r)
          for (int64_t d = 0; d < n_targets; ++d) target_idx[r * n_targets + d] = 2 * (target_n(r, d) - n_min_arr[r]);

      } else {
        NDA_RUNTIME_ERROR << "nfft_buf_t: only type3 and direct supported with target frequencies\n";
      }
    }

    /**
     * Convenience constructor for Rank=1: accepts vector of matsubara_freq directly
     *
     * Same as the non-uniform target constructor, but for single-frequency objects (Rank=1).
     * Accepts a flat vector of matsubara_freq instead of vector<array<matsubara_freq, 1>>.
     */
    nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<mesh::matsubara_freq> const &target_mf_, int buf_size_, nfft_type_t type,
               double tol_ = 1e-15)
      requires(Rank == 1)
       : nfft_buf_t(std::move(fiw_vec_), to_array_vector(target_mf_), buf_size_, type, tol_) {}

    ~nfft_buf_t() {
      if (buf_counter != 0) std::cout << " WARNING: Points in NFFT Buffer lost \n";
      // plan automatically destroyed by unique_ptr
    }

    // nfft_buffer is move-only (unique_ptr member)
    nfft_buf_t(nfft_buf_t const &)            = delete;
    nfft_buf_t &operator=(nfft_buf_t const &) = delete;
    nfft_buf_t(nfft_buf_t &&)                 = default;
    nfft_buf_t &operator=(nfft_buf_t &&rhs) noexcept {
      // Custom move assignment: array_view::operator= does deep copy,
      // but we need to rebind views to point to the same underlying data.
      // Leverage working move constructor via destroy + placement new.
      if (this != &rhs) {
        std::destroy_at(this);
        std::construct_at(this, std::move(rhs));
      }
      return *this;
    }

    /// Rebind nfft buffer to new accumulation container of same shape
    void rebind(array_view<dcomplex, Rank> new_fiw_arr) {
      flush();
      TRIQS_ASSERT((new_fiw_arr.shape() == fiw_arr.shape() or fiw_arr.empty())
                   and " Nfft Buffer: Rebind to array of different shape not allowed ");
      fiw_arr.rebind(new_fiw_arr);
    }

    /// Insert tau-vector {tau_1, tau_2, ... } \in [0,\beta)^Rank and corresponding f(tau) into the NFFT buffer
    void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {

      // Check if buffer has been properly initialized
      if (x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      if (nfft_type == nfft_type_t::type1) {
        // Type 1: normalize tau to [-PI, PI) and apply phase correction
        double tau_sum = 0.0;
        for (int r = 0; r < Rank; ++r) {
          x_arr(r, buf_counter) = 2 * M_PI * (tau_arr[r] / beta - 0.5); // \in [-PI, PI)
          tau_sum += tau_arr[r];
        }
        fx_arr[buf_counter] = std::exp(dcomplex(0, M_PI * tau_sum / beta)) * ftau;
      } else {
        // Type 3 / direct: raw tau values, no normalization or phase correction needed
        for (int r = 0; r < Rank; ++r) x_arr(r, buf_counter) = tau_arr[r];
        fx_arr[buf_counter] = ftau;
      }

      ++buf_counter;

      // If buffer is full, perform transform
      if (is_full()) {
        do_nfft();
        buf_counter = 0;
      }
    }

    /// Flush contents of the nfft buffer
    void flush() {

      // Check if buffer has been properly initialized
      if (x_arr.empty()) NDA_RUNTIME_ERROR << " Using a default-constructed NFFT Buffer is not allowed\n";

      // Don't do anything if buffer is empty
      if (is_empty()) return;

      // Execute the transform
      do_nfft();
      buf_counter = 0;
    }

    private:
    // Transform type
    nfft_type_t nfft_type = nfft_type_t::type1;

    // Type 1: output array in matsubara frequencies (uniform grid)
    nda::array_view<dcomplex, Rank> fiw_arr;

    // 1D output vector at target frequencies (type 3 and direct)
    nda::array_view<dcomplex, 1> fiw_vec;

    // Dimensions of the output array (type 1 only)
    std::array<int64_t, Rank> niws{};

    // Finufft plan (RAII-managed)
    finufft_plan_ptr plan;

    // Number of tau points for the nfft
    int buf_size = 0;

    // Inverse temperature (type 1 and direct)
    double beta = 0;

    // Counter for elements currently in the buffer
    int buf_counter = 0;

    // Common factor in container assignment (type 1 only)
    int common_factor = 1;

    // Number of target frequencies (type 3 and direct)
    int64_t n_targets = 0;

    // FINUFFT options struct
    finufft_opts opts{};

    // Array containing x values for the NFFT transform
    nda::array<double, 2> x_arr;

    // Array containing f(x) values for the NFFT transform
    nda::vector<dcomplex> fx_arr;

    // Array containing the NFFT output h(k) (type 1)
    nda::array<dcomplex, Rank> fk_arr;

    // Target frequencies, shape (Rank, n_targets) (type 3 only)
    nda::array<double, 2> s_arr;

    // Type 3 NFFT output buffer
    nda::vector<dcomplex> fk_vec;

    // Integer Matsubara indices per target, shape (Rank, n_targets) (direct only)
    nda::array<long, 2> target_n;

    // Per-dimension min index and range for power table (direct only)
    std::array<long, Rank> n_min_arr{};
    std::array<long, Rank> n_range_arr{};

    // Preallocated power table for direct mode (avoids repeated allocation)
    // Shape: (buf_size, n_range_arr[r]) per rank - stores powers for all buffer elements
    mutable std::array<nda::array<dcomplex, 2>, Rank> pow_tbl;

    // Precomputed offset indices for direct mode, doubled for AoS gather:
    // target_idx[r * n_targets + d] = 2 * (target_n(r, d) - n_min_arr[r])
    std::vector<long> target_idx;

    // Tolerance for the transformation
    double tol = 1e-15;

    // Function to check whether buffer is filled
    bool is_full() const { return buf_counter >= buf_size; }

    // Function to check whether buffer is empty
    bool is_empty() const { return buf_counter == 0; }

    // Perform NFFT transform and accumulate
    void do_nfft() {
      if (nfft_type == nfft_type_t::type1)
        do_nfft_type1();
      else if (nfft_type == nfft_type_t::type3)
        do_nfft_type3();
      else
        do_direct();
    }

    // Set source points and optional target points on the FINUFFT plan, dispatching by Rank.
    // FINUFFT expects coordinates in reverse rank order (x=last dim, y=second-to-last, z=first).
    // For type 1, pass nullptr for tgt. For type 3, pass target coordinate array.
    void set_pts(nda::array<double, 2> *tgt = nullptr) {
      auto _ = nda::range::all;
      auto n_tgt = tgt ? n_targets : int64_t{0};
      auto t = [&](int r) -> double * { return tgt ? (*tgt)(r, _).data() : nullptr; };
      if constexpr (Rank == 1)
        check_finufft(finufft_setpts(plan.get(), buf_counter, x_arr(0, _).data(), nullptr, nullptr, n_tgt, t(0), nullptr, nullptr));
      else if constexpr (Rank == 2)
        check_finufft(finufft_setpts(plan.get(), buf_counter, x_arr(1, _).data(), x_arr(0, _).data(), nullptr, n_tgt, t(1), t(0), nullptr));
      else // Rank == 3
        check_finufft(finufft_setpts(plan.get(), buf_counter, x_arr(2, _).data(), x_arr(1, _).data(), x_arr(0, _).data(), n_tgt, t(2), t(1), t(0)));
    }

    // Type 1: non-uniform tau -> uniform grid
    void do_nfft_type1() {
      set_pts();
      check_finufft(finufft_execute(plan.get(), fx_arr.data(), fk_arr.data()));

      // Accumulate results in fiw_arr. Care to normalize results afterwards
      for (auto idx_tpl : fiw_arr.indices()) {
        auto idx_sum = std::apply([](auto... idx) { return (idx + ... + 0); }, idx_tpl);
        int factor   = common_factor * (idx_sum % 2 ? -1 : 1);
        std::apply(fiw_arr, idx_tpl) += std::apply(fk_arr, idx_tpl) * factor;
      }
    }

    // Type 3: non-uniform tau -> non-uniform target frequencies
    void do_nfft_type3() {
      set_pts(&s_arr);
      check_finufft(finufft_execute(plan.get(), fx_arr.data(), fk_vec.data()));

      fiw_vec += fk_vec;
    }

    // Direct DFT: exploits Matsubara structure exp(i*omega_n*tau) = z^(2n+1)
    void do_direct() {
      using cbatch                    = xsimd::batch<dcomplex>;
      using rbatch                    = xsimd::batch<double>;
      using ibatch                    = xsimd::batch<int64_t>;
      constexpr std::size_t simd_size = cbatch::size;

      double const pi_over_beta    = M_PI / beta;
      int64_t const n_targets_simd = n_targets - (n_targets % static_cast<int64_t>(simd_size));
      dcomplex *fiw_ptr            = fiw_vec.data();

      // Phase 1: Build power tables pow_tbl[r](j,i) = z^(2*(n_min+i)+1)
      for (int r = 0; r < Rank; ++r) {
        long const n_range      = n_range_arr[r];
        long const n_range_simd = n_range - (n_range % static_cast<long>(simd_size));
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
          alignas(cbatch::arch_type::alignment()) std::array<dcomplex, simd_size> mult_arr;
          dcomplex z2_pow = 1.0;
          for (std::size_t k = 0; k < simd_size; ++k) {
            mult_arr[k] = z2_pow;
            z2_pow *= z2;
          }
          cbatch const mult(cbatch::load_aligned(mult_arr.data()));
          cbatch const stride(z2_pow);
          cbatch zp_vec(zp);

          for (long i = 0; i < n_range_simd; i += simd_size) {
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
        poet::static_for<0, Rank>([&](auto r) {
          pow_ptr[r] = reinterpret_cast<double const *>(&pow_tbl[r](j, 0));
        });

        // SIMD loop with gather
        for (int64_t d = 0; d < n_targets_simd; d += simd_size) {
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

  // Deduction guides for nfft_buf_t

  // Type 1: deduce Rank from output array
  template <int Rank>
  nfft_buf_t(nda::array_view<dcomplex, Rank>, int, double, double) -> nfft_buf_t<Rank>;

  // Type 3/direct: deduce Rank from target frequency array
  template <std::size_t N>
  nfft_buf_t(nda::array_view<dcomplex, 1>, std::vector<std::array<mesh::matsubara_freq, N>>, int, nfft_type_t, double) -> nfft_buf_t<static_cast<int>(N)>;

  // Convenience: vector<matsubara_freq> implies Rank=1
  nfft_buf_t(nda::array_view<dcomplex, 1>, std::vector<mesh::matsubara_freq> const &, int, nfft_type_t, double) -> nfft_buf_t<1>;

} // namespace triqs::utility
