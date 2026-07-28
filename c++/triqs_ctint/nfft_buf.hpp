// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <nda/nda.hpp>
#include <array>
#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>
#include <triqs/mesh/matsubara_freq.hpp>
#include <triqs/utility/exceptions.hpp> // TRIQS_ASSERT

namespace triqs::utility {

  namespace detail {
    // Opaque handle for the FINUFFT plan, defined in nfft_buf.cpp so that using this header
    // requires no FINUFFT header. The stateless deleter is empty-base-optimised, so the
    // unique_ptr member below stays pointer-sized, and its out-of-line operator() sidesteps
    // default_delete's requirement that the pointee be a complete type.
    struct nfft_plan;
    struct nfft_plan_deleter {
      void operator()(nfft_plan *) const;
    };
  } // namespace detail

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
    nfft_buf_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, double tol_ = 1e-15);

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
     * @param tol_        FINUFFT tolerance for type3 (default 1e-13). Ignored for direct.
     *                    Note: at upsampfac=2 FINUFFT caps the kernel width, so type3 cannot
     *                    reach 1e-15 (would need ns=17 > 16); 1e-13 is the practical default.
     */
    nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<std::array<mesh::matsubara_freq, Rank>> target_mf_, int buf_size_, nfft_type_t type,
               double tol_ = 1e-13);

    /**
     * Convenience constructor for Rank=1: accepts vector of matsubara_freq directly
     *
     * Same as the non-uniform target constructor, but for single-frequency objects (Rank=1).
     * Accepts a flat vector of matsubara_freq instead of vector<array<matsubara_freq, 1>>.
     */
    nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<mesh::matsubara_freq> const &target_mf_, int buf_size_, nfft_type_t type,
               double tol_ = 1e-13)
      requires(Rank == 1)
       : nfft_buf_t(std::move(fiw_vec_), to_array_vector(target_mf_), buf_size_, type, tol_) {}

    ~nfft_buf_t();

    // nfft_buffer is move-only (unique_ptr member)
    nfft_buf_t(nfft_buf_t const &)            = delete;
    nfft_buf_t &operator=(nfft_buf_t const &) = delete;
    nfft_buf_t(nfft_buf_t &&)                 = default;
    nfft_buf_t &operator=(nfft_buf_t &&rhs) noexcept;

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
    std::unique_ptr<detail::nfft_plan, detail::nfft_plan_deleter> plan;

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

    // Binary exponentiation approach: store z^(2^k) for k=0,1,2,...
    int num_power2_levels = 0;

    // Preallocated power-of-2 table for Rank==1 bitwise direct kernel.
    // Shape: (num_power2_levels, buf_size).
    // pow2_tbl(k, j) = z^(2^k) for buffer element j.
    std::conditional_t<Rank == 1, nda::array<dcomplex, 2>, std::monostate> pow2_tbl;

    // Per-target list of power-of-two exponents for Rank==1 bitwise kernel.
    // target_pow2_bits[d] contains k such that |2*n_d+1| has bit k set.
    std::conditional_t<Rank == 1, std::vector<std::vector<int>>, std::monostate> target_pow2_bits;

    // Rank>1 prime-sum direct kernel data
    std::vector<int> primes;
    std::array<nda::array<dcomplex, 2>, Rank> prime_pow_tbl;
    std::array<std::vector<std::vector<int>>, Rank> target_prime_sums;

    // Tolerance for the transformation
    double tol = 1e-15;

    // Function to check whether buffer is filled
    bool is_full() const { return buf_counter >= buf_size; }

    // Function to check whether buffer is empty
    bool is_empty() const { return buf_counter == 0; }

    // For fermionic frequencies: omega_n = (2n+1) * pi / beta.
    // Direct kernels always work with the absolute odd exponent |2n+1|.
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

    static constexpr int max_prime_sum_terms = 8;
    static constexpr int prime_sum_precompute_size = 128;

    // Tunable number of SIMD accumulators for direct kernels.
    // `n_acc_bitwise` is used by the Rank-1 bitwise power-of-two kernel.
    // `n_acc_prime` is used by the Rank>1 prime-sum kernel.
    static constexpr int n_acc_bitwise = 4;
    static constexpr int n_acc_prime = 4;
    // for rank 2 large sizes 8 is better but I think finufft will be faster anyway at that point

    struct prime_sum_entry_t {
      std::array<int, max_prime_sum_terms> terms{};
      int size = 0;
    };

    static constexpr prime_sum_entry_t express_as_prime_sum_ct(long n) {
      prime_sum_entry_t out{};
      while (n > 0 && out.size < max_prime_sum_terms) {
        if (n == 1) {
          out.terms[out.size++] = 1;
          break;
        }
        if (n == 2 || n == 3) {
          out.terms[out.size++] = static_cast<int>(n);
          break;
        }
        if (n == 4) {
          out.terms[out.size++] = 2;
          out.terms[out.size++] = 2;
          break;
        }

        long p = n;
        while (p > 1 && !is_prime(p)) --p;
        out.terms[out.size++] = static_cast<int>(p);
        n -= p;
      }
      return out;
    }

    static constexpr auto precomputed_prime_sums = [] {
      std::array<prime_sum_entry_t, prime_sum_precompute_size> table{};
      for (int n = 0; n < prime_sum_precompute_size; ++n) table[n] = express_as_prime_sum_ct(n);
      return table;
    }();

    static constexpr bool check_precomputed_prime_sums() {
      for (int n = 0; n < prime_sum_precompute_size; ++n) {
        long sum = 0;
        for (int i = 0; i < precomputed_prime_sums[n].size; ++i) sum += precomputed_prime_sums[n].terms[i];
        if (sum != n) return false;
      }
      return true;
    }

    static_assert(check_precomputed_prime_sums(), "prime-sum precompute table is invalid");

    // Helper: express n as sum of primes with repetition: n = p1 + p2 + ... + pk.
    static std::vector<int> express_as_prime_sum(long n) {
      if (n < 1) return {};
      auto entry = express_as_prime_sum_ct(n);
      return {entry.terms.begin(), entry.terms.begin() + entry.size};
    }

    // Accumulates n_acc independent targets at once: SIMD over buffer elements, instruction-level
    // parallelism across the targets. Defined in nfft_buf.cpp, where the SIMD types live.
    template <int n_acc, typename SimdPowFunc, typename ScalarPowFunc>
    [[gnu::always_inline]] inline void accumulate_targets_ilp(int64_t n_targets_total, int64_t buf_counter_simd, dcomplex *fiw_ptr,
                                                             SimdPowFunc &&compute_simd_pow, ScalarPowFunc &&compute_scalar_pow);

    // Perform NFFT transform and accumulate
    void do_nfft();

    // Set source points and optional target points on the FINUFFT plan, dispatching by Rank.
    // FINUFFT expects coordinates in reverse rank order (x=last dim, y=second-to-last, z=first).
    // For type 1, pass nullptr for tgt. For type 3, pass target coordinate array.
    void set_pts(nda::array<double, 2> *tgt = nullptr);

    // Type 1: non-uniform tau -> uniform grid
    void do_nfft_type1();

    // Type 3: non-uniform tau -> non-uniform target frequencies
    void do_nfft_type3();

    // Rank-1 direct NUDFT via bitwise power-of-two decomposition.
    // Constrained rather than asserted: the explicit instantiations in nfft_buf.cpp
    // instantiate every unconstrained member, and this body is valid for Rank 1 only.
    void do_direct_bitwise()
      requires(Rank == 1);

    // Rank>1 direct NUDFT via prime-sum decomposition
    void do_direct_prime();

    // Direct DFT dispatcher: bitwise decomposition for Rank 1, prime-sum for Rank>1
    void do_direct();
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
