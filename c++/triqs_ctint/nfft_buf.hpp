// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <nda/nda.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <vector>

#include <triqs/mesh/matsubara_freq.hpp>

#include "finufft.h"

namespace triqs::utility {

  inline void check_finufft(int err) {
    if (err > 0) NDA_RUNTIME_ERROR << "Error in FINUFFT: " << err << "\n";
  }

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

    /// Type 1 constructor: non-uniform tau -> uniform Matsubara grid
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
      check_finufft(finufft_makeplan(/*type =*/1, Rank, Ns.data(), /*iflag=*/1, /*n_transf =*/1, tol, &plan, &opts));
    }

    /// Non-uniform target constructor: type3 (FINUFFT) or direct DFT
    /// target_mf: vector of target frequency points (each point is an array of Rank matsubara_freq)
    /// type: nfft_type_t::type3 or nfft_type_t::direct
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
        check_finufft(finufft_makeplan(3, Rank, nullptr, /*iflag=*/1, /*n_transf=*/1, tol, &plan, &opts));

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
        // Preallocate power tables
        for (int r = 0; r < Rank; ++r) pow_tbl[r].resize(n_range_arr[r]);

        // Precompute offset indices: target_idx[r * n_targets + d] = target_n(r, d) - n_min_arr[r]
        target_idx.resize(Rank * n_targets);
        for (int r = 0; r < Rank; ++r)
          for (int64_t d = 0; d < n_targets; ++d) target_idx[r * n_targets + d] = target_n(r, d) - n_min_arr[r];

      } else {
        NDA_RUNTIME_ERROR << "nfft_buf_t: only type3 and direct supported with target frequencies\n";
      }
    }

    /// Convenience constructor for Rank=1: accepts vector of matsubara_freq directly
    nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<mesh::matsubara_freq> const &target_mf_, int buf_size_, nfft_type_t type,
               double tol_ = 1e-15)
      requires(Rank == 1)
       : nfft_buf_t(std::move(fiw_vec_), to_array_vector(target_mf_), buf_size_, type, tol_) {}

    ~nfft_buf_t() {
      if (buf_counter != 0) std::cout << " WARNING: Points in NFFT Buffer lost \n";
      if (plan) finufft_destroy(plan);
    }

    // nfft_buffer needs to be uncopyable, because nfft_plan contains raw pointers
    nfft_buf_t(nfft_buf_t const &)            = delete;
    nfft_buf_t(nfft_buf_t &&)                 = default;
    nfft_buf_t &operator=(nfft_buf_t const &) = delete;
    nfft_buf_t &operator=(nfft_buf_t &&rhs) noexcept {
      nfft_type = rhs.nfft_type;
      fiw_arr.rebind(rhs.fiw_arr);
      fiw_vec.rebind(rhs.fiw_vec);
      niws = rhs.niws;
      std::swap(plan, rhs.plan);
      buf_size      = rhs.buf_size;
      beta          = rhs.beta;
      buf_counter   = rhs.buf_counter;
      common_factor = rhs.common_factor;
      n_targets     = rhs.n_targets;
      opts          = std::move(rhs.opts);
      x_arr         = std::move(rhs.x_arr);
      fx_arr        = std::move(rhs.fx_arr);
      fk_arr        = std::move(rhs.fk_arr);
      s_arr         = std::move(rhs.s_arr);
      fk_vec        = std::move(rhs.fk_vec);
      target_n      = std::move(rhs.target_n);
      n_min_arr     = rhs.n_min_arr;
      n_range_arr   = rhs.n_range_arr;
      pow_tbl       = std::move(rhs.pow_tbl);
      target_idx    = std::move(rhs.target_idx);
      tol           = rhs.tol;
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

    // Finufft plan
    finufft_plan plan{nullptr};

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
    mutable std::array<std::vector<dcomplex>, Rank> pow_tbl;

    // Precomputed offset indices for direct mode: target_idx[r * n_targets + d] = target_n(r, d) - n_min_arr[r]
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
        check_finufft(finufft_setpts(plan, buf_counter, x_arr(0, _).data(), nullptr, nullptr, n_tgt, t(0), nullptr, nullptr));
      else if constexpr (Rank == 2)
        check_finufft(finufft_setpts(plan, buf_counter, x_arr(1, _).data(), x_arr(0, _).data(), nullptr, n_tgt, t(1), t(0), nullptr));
      else // Rank == 3
        check_finufft(finufft_setpts(plan, buf_counter, x_arr(2, _).data(), x_arr(1, _).data(), x_arr(0, _).data(), n_tgt, t(2), t(1), t(0)));
    }

    // Type 1: non-uniform tau -> uniform grid
    void do_nfft_type1() {
      set_pts();
      check_finufft(finufft_execute(plan, fx_arr.data(), fk_arr.data()));

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
      check_finufft(finufft_execute(plan, fx_arr.data(), fk_vec.data()));

      fiw_vec += fk_vec;
    }

    // Direct DFT: exploits Matsubara structure exp(i*omega_n*tau) = z^(2n+1)
    void do_direct() {

      // Get raw pointers for inner loop (avoid bounds checking)
      std::array<dcomplex *, Rank> pow_ptr;
      for (int r = 0; r < Rank; ++r) pow_ptr[r] = pow_tbl[r].data();
      dcomplex *fiw_ptr = fiw_vec.data();
      long const *idx_ptr = target_idx.data();

      // Precompute invariants
      double const pi_over_beta = M_PI / beta;
      std::array<long, Rank> abs_n_min;
      std::array<bool, Rank> n_min_positive;
      for (int r = 0; r < Rank; ++r) {
        abs_n_min[r] = std::abs(n_min_arr[r]);
        n_min_positive[r] = n_min_arr[r] >= 0;
      }

      for (int j = 0; j < buf_counter; ++j) {
        // Build power table z^(2n+1) for each dimension
        for (int r = 0; r < Rank; ++r) {
          double theta = pi_over_beta * x_arr(r, j);
          dcomplex z{std::cos(theta), std::sin(theta)};
          dcomplex z2 = z * z;

          // Compute z^(2*n_min+1) = z * z2^n_min, using conj(z2) = 1/z2 since |z|=1
          dcomplex zp  = z;
          dcomplex z2n = n_min_positive[r] ? z2 : std::conj(z2);
          for (long i = 0; i < abs_n_min[r]; ++i) zp *= z2n;

          pow_ptr[r][0] = zp;
          for (long i = 1; i < n_range_arr[r]; ++i) {
            zp *= z2;
            pow_ptr[r][i] = zp;
          }
        }

        // Accumulate using precomputed powers and raw pointers (unrolled by rank)
        dcomplex fj = fx_arr[j];
        if constexpr (Rank == 1) {
          for (int64_t d = 0; d < n_targets; ++d) fiw_ptr[d] += fj * pow_ptr[0][idx_ptr[d]];
        } else if constexpr (Rank == 2) {
          long const *idx0 = idx_ptr;
          long const *idx1 = idx_ptr + n_targets;
          for (int64_t d = 0; d < n_targets; ++d) fiw_ptr[d] += fj * pow_ptr[0][idx0[d]] * pow_ptr[1][idx1[d]];
        } else { // Rank == 3
          long const *idx0 = idx_ptr;
          long const *idx1 = idx_ptr + n_targets;
          long const *idx2 = idx_ptr + 2 * n_targets;
          for (int64_t d = 0; d < n_targets; ++d) fiw_ptr[d] += fj * pow_ptr[0][idx0[d]] * pow_ptr[1][idx1[d]] * pow_ptr[2][idx2[d]];
        }
      }
    }
  };
} // namespace triqs::utility
