// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

// Everything that needs FINUFFT or xsimd lives here rather than in nfft_buf.hpp, so that neither
// library appears in the installed interface of this project. nfft_buf_t is explicitly
// instantiated for all Ranks it supports at the bottom of this file.

#include "./nfft_buf.hpp"

#include "finufft.h"
#include <xsimd/xsimd.hpp>

#include <algorithm>
#include <bit>
#include <iostream>
#include <utility>

namespace triqs::utility {

  namespace {
    void check_finufft(int err) {
      if (err > 0) NDA_RUNTIME_ERROR << "Error in FINUFFT: " << err << "\n";
    }
  } // namespace

  namespace detail {

    // Compile-time unrolled loop over [0, N): invokes f with std::integral_constant<int, I> for each
    // I. The index is int, matching the int Rank / n_acc bounds and keeping target arithmetic signed.
    template <int N, typename F> constexpr void static_for(F &&f) {
      [&]<std::size_t... Is>(std::index_sequence<Is...>) { (f(std::integral_constant<int, Is>{}), ...); }(std::make_index_sequence<N>{});
    }

    // The opaque handle declared in nfft_buf.hpp: finufft_plan is itself a pointer to an
    // incomplete type, so it needs a definition of our own to hang the deleter off.
    struct nfft_plan {
      finufft_plan p = nullptr;
    };

    void nfft_plan_deleter::operator()(nfft_plan *ptr) const {
      if (ptr) {
        if (ptr->p) finufft_destroy(ptr->p);
        delete ptr;
      }
    }

  } // namespace detail

  // ═══════════════════════════════════════════════════════════════════════════
  // Exponent decomposition for the direct kernels
  //
  // None of this depends on Rank, so it lives here as file-local free functions
  // rather than as static members of nfft_buf_t. Used only when setting up the
  // per-target exponent lists in the non-uniform-target constructor.
  // ═══════════════════════════════════════════════════════════════════════════

  namespace {

    // For fermionic frequencies: omega_n = (2n+1) * pi / beta.
    // Direct kernels always work with the absolute odd exponent |2n+1|.
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

    constexpr int max_prime_sum_terms       = 8;
    constexpr int prime_sum_precompute_size = 128;

    // Tunable number of SIMD accumulators for direct kernels.
    // `n_acc_bitwise` is used by the Rank-1 bitwise power-of-two kernel.
    // `n_acc_prime` is used by the Rank>1 prime-sum kernel.
    constexpr int n_acc_bitwise = 4;
    constexpr int n_acc_prime   = 4;
    // for rank 2 large sizes 8 is better but I think finufft will be faster anyway at that point

    struct prime_sum_entry_t {
      std::array<int, max_prime_sum_terms> terms{};
      int size = 0;
    };

    constexpr prime_sum_entry_t express_as_prime_sum_ct(long n) {
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

    // Compile-time self-check of express_as_prime_sum_ct: every decomposition must sum back to n.
    // The table exists only to drive the static_assert below; the runtime helper recomputes.
    constexpr auto precomputed_prime_sums = [] {
      std::array<prime_sum_entry_t, prime_sum_precompute_size> table{};
      for (int n = 0; n < prime_sum_precompute_size; ++n) table[n] = express_as_prime_sum_ct(n);
      return table;
    }();

    constexpr bool check_precomputed_prime_sums() {
      for (int n = 0; n < prime_sum_precompute_size; ++n) {
        long sum = 0;
        for (int i = 0; i < precomputed_prime_sums[n].size; ++i) sum += precomputed_prime_sums[n].terms[i];
        if (sum != n) return false;
      }
      return true;
    }

    static_assert(check_precomputed_prime_sums(), "prime-sum precompute table is invalid");

    // Helper: express n as sum of primes with repetition: n = p1 + p2 + ... + pk.
    std::vector<int> express_as_prime_sum(long n) {
      if (n < 1) return {};
      auto entry = express_as_prime_sum_ct(n);
      return {entry.terms.begin(), entry.terms.begin() + entry.size};
    }

  } // namespace

  // ═══════════════════════════════════════════════════════════════════════════
  // Construction and destruction
  // ═══════════════════════════════════════════════════════════════════════════

  template <int Rank>
  nfft_buf_t<Rank>::nfft_buf_t(nda::array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, double tol_)
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

    // Init nfft_plan. finufft_makeplan copies what it needs from opts, so it is a local.
    finufft_opts opts{};
    finufft_default_opts(&opts); // set default opts (must start with this)
    opts.nthreads = 1;           // enforce single-thread
    auto Ns = std::vector(niws.rbegin(), niws.rend()); // Reverse order for FINUFFT
    finufft_plan raw_plan = nullptr;
    check_finufft(finufft_makeplan(/*type =*/1, Rank, Ns.data(), /*iflag=*/1, /*n_transf =*/1, tol, &raw_plan, &opts));
    plan.reset(new detail::nfft_plan{raw_plan});
  }

  template <int Rank>
  nfft_buf_t<Rank>::nfft_buf_t(nda::array_view<dcomplex, 1> fiw_vec_, std::vector<std::array<mesh::matsubara_freq, Rank>> target_mf_, int buf_size_,
                               nfft_type_t type, double tol_)
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
      finufft_opts opts{};
      finufft_default_opts(&opts);
      opts.nthreads = 1;
      finufft_plan raw_plan = nullptr;
      check_finufft(finufft_makeplan(3, Rank, nullptr, /*iflag=*/1, /*n_transf=*/1, tol, &raw_plan, &opts));
      plan.reset(new detail::nfft_plan{raw_plan});

    } else if (type == nfft_type_t::direct) {
      // Extract integer indices and beta from matsubara_freq
      beta = target_mf_[0][0].beta;
      target_n.resize(Rank, n_targets);
      for (int r = 0; r < Rank; ++r)
        for (int64_t d = 0; d < n_targets; ++d) target_n(r, d) = target_mf_[d][r].n;

      if constexpr (Rank > 1) {
        // Rank>1: use prime-sum direct kernel.
        std::vector<int> all_primes;
        for (int r = 0; r < Rank; ++r) {
          target_prime_sums[r].resize(n_targets);
          for (int64_t d = 0; d < n_targets; ++d) {
            auto exponent           = odd_exponent_abs(target_n(r, d));
            auto prime_list         = express_as_prime_sum(static_cast<long>(exponent));
            target_prime_sums[r][d] = prime_list;
            for (int prime : prime_list) all_primes.push_back(prime);
          }
        }

        std::sort(all_primes.begin(), all_primes.end());
        all_primes.erase(std::unique(all_primes.begin(), all_primes.end()), all_primes.end());
        primes = std::move(all_primes);

        for (int r = 0; r < Rank; ++r) {
          for (int64_t d = 0; d < n_targets; ++d) {
            for (int &prime : target_prime_sums[r][d]) {
              prime = static_cast<int>(std::find(primes.begin(), primes.end(), prime) - primes.begin());
            }
          }
        }

        for (int r = 0; r < Rank; ++r) { prime_pow_tbl[r].resize(primes.size(), buf_size_); }
      } else {
        // Rank-1: use bitwise power-of-two direct kernel.
        unsigned long max_exponent = 0;
        target_pow2_bits.resize(n_targets);
        for (int64_t d = 0; d < n_targets; ++d) {
          unsigned long exponent = odd_exponent_abs(target_n(0, d));
          max_exponent           = std::max(max_exponent, exponent);

          std::vector<int> bits;
          for (int k = 0; exponent > 0; ++k, exponent >>= 1) {
            if (exponent & 1ul) bits.push_back(k);
          }
          target_pow2_bits[d] = std::move(bits);
        }

        num_power2_levels = std::max(1, static_cast<int>(std::bit_width(max_exponent)));
        pow2_tbl.resize(num_power2_levels, buf_size_);
      }

    } else {
      NDA_RUNTIME_ERROR << "nfft_buf_t: only type3 and direct supported with target frequencies\n";
    }
  }

  template <int Rank> nfft_buf_t<Rank>::~nfft_buf_t() {
    if (buf_counter != 0) std::cout << " WARNING: Points in NFFT Buffer lost \n";
    // plan automatically destroyed by unique_ptr
  }

  template <int Rank> nfft_buf_t<Rank> &nfft_buf_t<Rank>::operator=(nfft_buf_t &&rhs) noexcept {
    // Custom move assignment: array_view::operator= does deep copy,
    // but we need to rebind views to point to the same underlying data.
    // Leverage working move constructor via destroy + placement new.
    if (this != &rhs) {
      std::destroy_at(this);
      std::construct_at(this, std::move(rhs));
    }
    return *this;
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // Unified Target Accumulation Template with Instruction-Level Parallelism
  // ═══════════════════════════════════════════════════════════════════════════
  //
  // This template implements a sophisticated two-level parallelism strategy to maximize CPU throughput:
  //
  // 1. **SIMD (Single Instruction Multiple Data) Parallelism:**
  //    - Vectorizes across buffer elements (tau points)
  //    - Process 2-4 complex numbers simultaneously per instruction (hardware dependent)
  //    - Uses xsimd library for portable SIMD abstractions
  //
  // 2. **ILP (Instruction-Level Parallelism):**
  //    - Process n_acc independent targets simultaneously
  //    - Each target maintains its own accumulator chain
  //    - Breaks data dependencies, allowing CPU to execute multiple operations in parallel
  //    - Exploits superscalar execution and out-of-order execution in modern CPUs
  //
  // **Why this matters:**
  // Without ILP, the CPU would wait for each FMA (fused multiply-add) to complete before
  // starting the next one due to data dependencies. With n_acc=4 independent chains,
  // the CPU can overlap execution, achieving ~4x higher throughput.
  //
  // **Parameters:**
  // - n_acc: Number of independent accumulators (typically 4-8)
  // - fx: Buffered f(tau_j) values, buf_counter entries
  // - buf_counter: Number of buffered elements
  // - buf_counter_simd: buf_counter floored to a multiple of the SIMD width
  // - compute_simd_pow: Lambda computing exp(iω*τ) for SIMD-aligned buffer indices
  // - compute_scalar_pow: Lambda computing exp(iω*τ) for scalar tail elements
  //
  // Nothing here depends on Rank or on any other nfft_buf_t member, so it is a file-local free
  // function rather than a member template, keeping the declaration out of the installed header.
  //
  namespace {

    template <int n_acc, typename SimdPowFunc, typename ScalarPowFunc>
    [[gnu::always_inline]] inline void accumulate_targets_ilp(dcomplex const *fx, int buf_counter, int64_t buf_counter_simd, int64_t n_targets_total,
                                                              dcomplex *fiw_ptr, SimdPowFunc &&compute_simd_pow, ScalarPowFunc &&compute_scalar_pow) {
      using cbatch                    = xsimd::batch<dcomplex>; // SIMD type for complex numbers
      constexpr std::size_t simd_size = cbatch::size;           // Typically 2-4 depending on CPU

      // ─────────────────────────────────────────────────────────────────────────
      // Helper: Single-target accumulation (used for remainder targets)
      // ─────────────────────────────────────────────────────────────────────────
      // Computes: fiw[d] += Σ_j f(tau_j) * exp(iω_d*tau_j)
      auto accumulate_one = [&](int64_t d) {
        // SIMD loop: process buffer in chunks of simd_size
        cbatch sum_vec(dcomplex{0, 0}); // Initialize SIMD accumulator to zero
        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          cbatch fj  = cbatch::load_unaligned(fx + j); // Load f(tau_j) [simd_size elements]
          cbatch pow = compute_simd_pow(d, j);         // Compute exp(iω_d*tau_j) [vectorized]
          sum_vec    = xsimd::fma(fj, pow, sum_vec);   // Fused multiply-add: sum += fj * pow
        }
        // Reduce SIMD vector to scalar by summing all lanes
        dcomplex sum = xsimd::reduce_add(sum_vec);

        // Scalar tail: handle remaining buffer elements that don't fit in SIMD
        for (int j = buf_counter_simd; j < buf_counter; ++j) {
          dcomplex pow = compute_scalar_pow(d, j);
          sum += fx[j] * pow;
        }

        fiw_ptr[d] += sum; // Accumulate into output
      };

      // ─────────────────────────────────────────────────────────────────────────
      // Main ILP Loop: Process n_acc targets simultaneously
      // ─────────────────────────────────────────────────────────────────────────
      // Round down to nearest multiple of n_acc
      int64_t const n_targets_main = (n_targets_total / n_acc) * n_acc;
      int64_t d                    = 0;

      for (; d < n_targets_main; d += n_acc) {
        // Initialize n_acc independent SIMD accumulators (one per target in this batch)
        std::array<cbatch, n_acc> sum_vecs;
        detail::static_for<n_acc>([&](const auto acc_idx) { sum_vecs[acc_idx] = cbatch(dcomplex{0, 0}); });

        // ═══ SIMD Loop: Vectorize over buffer elements ═══
        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          // Load f(tau_j) once - shared across all n_acc targets
          cbatch fj = cbatch::load_unaligned(fx + j);

          // Unroll over n_acc accumulators at compile time
          // This creates n_acc independent FMA dependency chains, enabling ILP
          // The CPU can execute these in parallel pipelines
          detail::static_for<n_acc>([&](const auto acc_idx) {
            cbatch pow        = compute_simd_pow(d + acc_idx, j);       // exp(iω_{d+k}*tau_j)
            sum_vecs[acc_idx] = xsimd::fma(fj, pow, sum_vecs[acc_idx]); // Independent accumulation
          });
        }

        // ═══ Reduce Phase: SIMD vectors → scalars ═══
        std::array<dcomplex, n_acc> sums;
        detail::static_for<n_acc>([&](const auto acc_idx) { sums[acc_idx] = xsimd::reduce_add(sum_vecs[acc_idx]); });

        // ═══ Scalar Tail: Process remaining non-SIMD-aligned elements ═══
        for (int j = buf_counter_simd; j < buf_counter; ++j) {
          dcomplex fj = fx[j];
          detail::static_for<n_acc>([&](const auto acc_idx) {
            dcomplex pow = compute_scalar_pow(d + acc_idx, j);
            sums[acc_idx] += fj * pow;
          });
        }

        // ═══ Write-back Phase: Store results to output array ═══
        detail::static_for<n_acc>([&](const auto acc_idx) { fiw_ptr[d + acc_idx] += sums[acc_idx]; });
      }

      // ─────────────────────────────────────────────────────────────────────────
      // Remainder Loop: Handle final targets when n_targets % n_acc != 0
      // ─────────────────────────────────────────────────────────────────────────
      for (; d < n_targets_total; ++d) accumulate_one(d);
    }

  } // namespace

  // ═══════════════════════════════════════════════════════════════════════════
  // FINUFFT-backed transforms
  // ═══════════════════════════════════════════════════════════════════════════

  template <int Rank> void nfft_buf_t<Rank>::set_pts(nda::array<double, 2> *tgt) {
    auto _ = nda::range::all;
    auto n_tgt = tgt ? n_targets : int64_t{0};
    auto t = [&](int r) -> double * { return tgt ? (*tgt)(r, _).data() : nullptr; };
    if constexpr (Rank == 1)
      check_finufft(finufft_setpts(plan->p, buf_counter, x_arr(0, _).data(), nullptr, nullptr, n_tgt, t(0), nullptr, nullptr));
    else if constexpr (Rank == 2)
      check_finufft(finufft_setpts(plan->p, buf_counter, x_arr(1, _).data(), x_arr(0, _).data(), nullptr, n_tgt, t(1), t(0), nullptr));
    else // Rank == 3
      check_finufft(finufft_setpts(plan->p, buf_counter, x_arr(2, _).data(), x_arr(1, _).data(), x_arr(0, _).data(), n_tgt, t(2), t(1), t(0)));
  }

  template <int Rank> void nfft_buf_t<Rank>::do_nfft_type1() {
    set_pts();
    check_finufft(finufft_execute(plan->p, fx_arr.data(), fk_arr.data()));

    // Accumulate results in fiw_arr. Care to normalize results afterwards
    for (auto idx_tpl : fiw_arr.indices()) {
      auto idx_sum = std::apply([](auto... idx) { return (idx + ... + 0); }, idx_tpl);
      int factor   = common_factor * (idx_sum % 2 ? -1 : 1);
      std::apply(fiw_arr, idx_tpl) += std::apply(fk_arr, idx_tpl) * factor;
    }
  }

  template <int Rank> void nfft_buf_t<Rank>::do_nfft_type3() {
    set_pts(&s_arr);
    check_finufft(finufft_execute(plan->p, fx_arr.data(), fk_vec.data()));

    fiw_vec += fk_vec;
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // Rank-1 Direct NUDFT: Bitwise Power-of-Two Decomposition
  // ═══════════════════════════════════════════════════════════════════════════
  //
  // **Algorithm Overview:**
  // Goal: Compute exp(iω_n*τ) for fermionic Matsubara frequencies ω_n = (2n+1)π/β
  //
  // **Key Mathematical Insight:**
  //   exp(iω_n*τ) = exp(i*(2n+1)*π*τ/β)
  //               = [exp(iπτ/β)]^(2n+1)
  //               = z^(2n+1)
  // where z := exp(iπτ/β) is the "base phase factor"
  //
  // **Efficient Exponentiation via Binary Decomposition:**
  // Instead of computing z^m naively (O(m) multiplications), we use binary exponentiation:
  //
  // 1. Express m = |2n+1| in binary: m = Σ b_k * 2^k  (where b_k ∈ {0,1} are the bits)
  //    Example: m=13 = 8+4+1 = 2³ + 2² + 2⁰
  //
  // 2. Precompute powers-of-two: z, z², z⁴, z⁸, z¹⁶, ... via repeated squaring
  //    This takes O(log m) operations
  //
  // 3. Multiply only the powers corresponding to set bits:
  //    z^m = z^(2^k₁) * z^(2^k₂) * ... where k₁, k₂, ... are the bit positions
  //    Example: z¹³ = z⁸ * z⁴ * z¹
  //
  // **Complexity:** O(log m) multiplications instead of O(m)
  // For typical Matsubara indices (|2n+1| ~ 1-1000), this is 10-100x faster!
  //
  // **Why this is optimal for Rank=1:**
  // Binary representation is minimal - every integer has exactly one binary form.
  // For Rank>1, we use prime-sum decomposition instead (better power sharing).
  //
  template <int Rank>
  void nfft_buf_t<Rank>::do_direct_bitwise()
    requires(Rank == 1)
  {
    using cbatch                    = xsimd::batch<dcomplex>;
    constexpr std::size_t simd_size = cbatch::size;

    double const pi_over_beta      = M_PI / beta;
    int64_t const buf_counter_simd = buf_counter & -simd_size;  // Floor to SIMD alignment
    dcomplex *fiw_ptr              = fiw_vec.data();            // Output pointer

    // ═══════════════════════════════════════════════════════════════════════
    // Phase 1: Build Power-of-Two Table via Repeated Squaring
    // ═══════════════════════════════════════════════════════════════════════
    // Compute pow2_tbl(k, j) = z_j^(2^k) for all buffer elements j
    // where z_j = exp(iπτ_j/β)

    // ─── Step 1a: Compute z = exp(iπτ/β) for all buffer elements ───
    // SIMD path: Process simd_size elements at once
    for (int j = 0; j < buf_counter_simd; j += simd_size) {
      using rbatch = xsimd::batch<double>;
      // Compute angle: θ = πτ/β
      rbatch theta_vec = rbatch::load_unaligned(&x_arr(0, j)) * pi_over_beta;
      // Vectorized sincos: compute sin(θ) and cos(θ) simultaneously (hardware optimized)
      auto [sin_vec, cos_vec] = xsimd::sincos(theta_vec);
      // Build complex exponential: z = cos(θ) + i*sin(θ) = exp(iθ)
      cbatch z_vec(cos_vec, sin_vec);
      // Store z^(2^0) = z in level k=0
      z_vec.store_unaligned(&pow2_tbl(0, j));
    }
    // Scalar tail: handle remaining elements that don't fit in SIMD
    for (int j = buf_counter_simd; j < buf_counter; ++j) {
      double const theta = pi_over_beta * x_arr(0, j);
      pow2_tbl(0, j)     = dcomplex{std::cos(theta), std::sin(theta)};
    }

    // ─── Step 1b: Repeated squaring to build higher levels ───
    // For each level k: z^(2^k) = [z^(2^(k-1))]²
    // Example: z² = z*z, z⁴ = z²*z², z⁸ = z⁴*z⁴, ...
    for (int k = 1; k < num_power2_levels; ++k) {
      // SIMD path: vectorized squaring
      for (int j = 0; j < buf_counter_simd; j += simd_size) {
        cbatch prev = cbatch::load_unaligned(&pow2_tbl(k - 1, j));  // Load z^(2^(k-1))
        cbatch curr = prev * prev;                                   // Square it
        curr.store_unaligned(&pow2_tbl(k, j));                      // Store z^(2^k)
      }
      // Scalar tail
      for (int j = buf_counter_simd; j < buf_counter; ++j) {
        dcomplex prev = pow2_tbl(k - 1, j);
        pow2_tbl(k, j) = prev * prev;
      }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Phase 2: Accumulate Targets by Multiplying Powers for Set Bits
    // ═══════════════════════════════════════════════════════════════════════
    // For each target d with exponent m = |2n_d+1|, we precomputed the bit
    // positions in target_pow2_bits[d]. Now we multiply those powers together.

    // ─── SIMD Power Computation Lambda ───
    auto compute_simd_pow = [&](int64_t d, int j) -> cbatch {
      auto const &bits  = target_pow2_bits[d];  // List of set bit positions in |2n_d+1|
      bool const is_neg = target_n(0, d) < 0;   // Is this a negative frequency?

      // Start with identity: z^0 = 1
      cbatch rank_pow(dcomplex{1.0, 0.0});

      // Multiply powers corresponding to each set bit
      // Example: if bits = [0, 2, 3], compute z^1 * z^4 * z^8 = z^13
      for (int k : bits) {
        rank_pow *= cbatch::load_unaligned(&pow2_tbl(k, j));
      }

      // For negative frequencies: exp(-iω*τ) = conj(exp(iω*τ))
      // This uses the identity: exp(-ix) = cos(x) - i*sin(x) = conj(exp(ix))
      return is_neg ? xsimd::conj(rank_pow) : rank_pow;
    };

    // ─── Scalar Power Computation Lambda (identical logic, non-vectorized) ───
    auto compute_scalar_pow = [&](int64_t d, int j) -> dcomplex {
      auto const &bits  = target_pow2_bits[d];
      bool const is_neg = target_n(0, d) < 0;
      dcomplex rank_pow{1.0, 0.0};
      for (int k : bits) rank_pow *= pow2_tbl(k, j);
      return is_neg ? std::conj(rank_pow) : rank_pow;
    };

    // ─── Call unified ILP accumulation template ───
    // This computes: fiw[d] += Σ_j f(tau_j) * exp(iω_d*tau_j) for all targets d
    accumulate_targets_ilp<n_acc_bitwise>(fx_arr.data(), buf_counter, buf_counter_simd, n_targets, fiw_ptr, compute_simd_pow, compute_scalar_pow);
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // Rank>1 Direct NUDFT: Prime-Sum Decomposition
  // ═══════════════════════════════════════════════════════════════════════════
  //
  // **Algorithm Overview:**
  // Goal: Compute exp(iω_{n1}*τ1 + iω_{n2}*τ2 + ...) for multi-dimensional frequencies
  //       = exp(iω_{n1}*τ1) * exp(iω_{n2}*τ2) * ...
  //       = z_1^{m1} * z_2^{m2} * ...
  // where z_r = exp(iπτ_r/β) and m_r = |2n_r+1|
  //
  // **Why Not Use Binary Exponentiation for Rank>1?**
  // For Rank>1, we need to compute many distinct exponents (m_r for each rank r and target d).
  // Binary decomposition would require storing 2^k powers for each distinct exponent,
  // leading to excessive memory usage and poor cache behavior.
  //
  // **Prime-Sum Decomposition Strategy:**
  // Every positive integer can be expressed as a sum of primes (Goldbach-style):
  //   m = p1 + p2 + ... + pk
  // Then: z^m = z^(p1+p2+...+pk) = z^p1 * z^p2 * ... * z^pk
  //
  // **Key Advantages:**
  // 1. **Power Sharing:** The set of unique primes needed across ALL targets and ranks
  //    is much smaller than the set of unique exponents. We compute each z^p once and reuse it.
  //
  // 2. **Small Exponents:** Primes are typically small (2, 3, 5, 7, 11, ...), making
  //    z^p fast to compute via binary exponentiation.
  //
  // 3. **Memory Efficiency:** Store O(#primes * Rank * buf_size) instead of
  //    O(#unique_exponents * Rank * buf_size). For many targets, #primes << #unique_exponents.
  //
  // **Example:**
  // Targets with m = 13, 15, 17 in some rank:
  //   13 = 13,        15 = 13+2,      17 = 17
  // Unique primes: {2, 13, 17}
  // We compute z^2, z^13, z^17 once, then:
  //   z^13 = z^13,    z^15 = z^13*z^2,  z^17 = z^17
  //
  template <int Rank> void nfft_buf_t<Rank>::do_direct_prime() {
    using cbatch                    = xsimd::batch<dcomplex>;
    constexpr std::size_t simd_size = cbatch::size;

    double const pi_over_beta      = M_PI / beta;
    int64_t const buf_counter_simd = buf_counter & -simd_size;
    dcomplex *fiw_ptr              = fiw_vec.data();
    int const num_primes           = static_cast<int>(primes.size());  // # of unique primes across all targets

    // ═══════════════════════════════════════════════════════════════════════
    // Phase 1: Compute Prime Powers for Each Rank
    // ═══════════════════════════════════════════════════════════════════════
    // For each rank r and each unique prime p, compute:
    //   prime_pow_tbl[r](p_idx, j) = z_r^p  where z_r = exp(iπτ_r/β)
    //
    // This is done once per rank, then reused for all targets.

    // Use compile-time loop over ranks (unrolled at compile time)
    detail::static_for<Rank>([&](const auto r) {
      // ─── Step 1: Compute base phase factors z_r for this rank ───
      std::vector<dcomplex> z_vals(buf_counter);
      for (int j = 0; j < buf_counter; ++j) {
        double const theta = pi_over_beta * x_arr(r, j);  // θ = πτ_r/β
        z_vals[j]          = dcomplex{std::cos(theta), std::sin(theta)};  // z_r = exp(iθ)
      }

      // ─── Step 2: For each unique prime, compute z_r^prime ───
      for (int p_idx = 0; p_idx < num_primes; ++p_idx) {
        int prime = primes[p_idx];

        // Special case: "prime" = 1 (treated as prime for algorithm simplicity)
        if (prime == 1) {
          // z^1 = z (no exponentiation needed)
          for (int j = 0; j < buf_counter; ++j) {
            prime_pow_tbl[r](p_idx, j) = z_vals[j];
          }
          continue;
        }

        // ═══ Binary Exponentiation to Compute z^prime ═══
        // This is O(log prime) instead of O(prime) for naive multiplication
        // Algorithm: Process bits of exponent from LSB to MSB
        //   result = 1
        //   base = z
        //   for each bit k in prime:
        //     if bit k is set: result *= base
        //     base = base^2  (square for next bit)

        // SIMD path: vectorized binary exponentiation
        for (int j = 0; j < buf_counter_simd; j += simd_size) {
          cbatch z_vec = cbatch::load_unaligned(&z_vals[j]);  // Load z
          cbatch zp_vec(dcomplex{1.0, 0.0});                   // Result accumulator (starts at 1)
          cbatch base_vec = z_vec;                             // Current power of base
          int exp         = prime;                             // Exponent to process

          // Binary exponentiation loop
          while (exp > 0) {
            if (exp & 1) zp_vec *= base_vec;  // If current bit is set, multiply into result
            base_vec *= base_vec;              // Square base for next bit position
            exp >>= 1;                         // Shift to next bit
          }
          zp_vec.store_unaligned(&prime_pow_tbl[r](p_idx, j));
        }

        // Scalar tail: same algorithm, non-vectorized
        for (int j = buf_counter_simd; j < buf_counter; ++j) {
          dcomplex zp{1.0, 0.0};
          dcomplex base = z_vals[j];
          int exp       = prime;
          while (exp > 0) {
            if (exp & 1) zp *= base;
            base *= base;
            exp >>= 1;
          }
          prime_pow_tbl[r](p_idx, j) = zp;
        }
      }
    });

    // ═══════════════════════════════════════════════════════════════════════
    // Phase 2: Accumulate Targets by Combining Prime Powers
    // ═══════════════════════════════════════════════════════════════════════
    // For each target d with multi-dimensional frequency (ω_{n1}, ω_{n2}, ...):
    //   exp(iω_{n1}*τ1 + ... + iω_{nR}*τR) = Π_r z_r^{m_r}
    // where m_r = |2n_r+1| and z_r = exp(iπτ_r/β)
    //
    // Each m_r is decomposed as sum of primes: m_r = p1 + p2 + ...
    // So: z_r^{m_r} = z_r^{p1} * z_r^{p2} * ...

    // ─── SIMD Power Computation Lambda ───
    auto compute_simd_pow = [&](int64_t d, int j) -> cbatch {
      cbatch pow_prod;  // Will accumulate product across all ranks

      // Loop over each dimension/rank (compile-time unroll)
      detail::static_for<Rank>([&](const auto r) {
        // Start with identity for this rank
        cbatch rank_pow(dcomplex{1.0, 0.0});

        // Get list of prime indices that sum to m_r = |2n_r+1|
        auto const &prime_indices = target_prime_sums[r][d];

        // Multiply prime powers: z_r^{m_r} = z_r^{p1} * z_r^{p2} * ...
        for (int prime_idx : prime_indices) {
          cbatch prime_pow = cbatch::load_unaligned(&prime_pow_tbl[r](prime_idx, j));
          rank_pow *= prime_pow;
        }

        // Handle negative frequencies via conjugation
        rank_pow = (target_n(r, d) < 0) ? xsimd::conj(rank_pow) : rank_pow;

        // Accumulate product across ranks: Π_r z_r^{m_r}
        pow_prod = (r == 0) ? rank_pow : pow_prod * rank_pow;
      });

      return pow_prod;
    };

    // ─── Scalar Power Computation Lambda (identical logic, non-vectorized) ───
    auto compute_scalar_pow = [&](int64_t d, int j) -> dcomplex {
      dcomplex pow_prod{1.0, 0.0};
      detail::static_for<Rank>([&](const auto r) {
        dcomplex rank_pow{1.0, 0.0};
        auto const &prime_indices = target_prime_sums[r][d];
        for (int prime_idx : prime_indices) {
          rank_pow *= prime_pow_tbl[r](prime_idx, j);
        }
        rank_pow  = (target_n(r, d) < 0) ? std::conj(rank_pow) : rank_pow;
        pow_prod *= rank_pow;
      });
      return pow_prod;
    };

    // ─── Call unified ILP accumulation template ───
    accumulate_targets_ilp<n_acc_prime>(fx_arr.data(), buf_counter, buf_counter_simd, n_targets, fiw_ptr, compute_simd_pow, compute_scalar_pow);
  }

  // ═══════════════════════════════════════════════════════════════════════════
  // Direct DFT Dispatcher: Choose Algorithm Based on Rank
  // ═══════════════════════════════════════════════════════════════════════════
  // Rank-1: Use bitwise power-of-two decomposition (optimal for single dimension)
  // Rank>1: Use prime-sum decomposition (better power sharing across dimensions)
  template <int Rank> void nfft_buf_t<Rank>::do_direct() {
    if constexpr (Rank == 1)
      do_direct_bitwise();
    else
      do_direct_prime();
  }

  // Perform NFFT transform and accumulate
  template <int Rank> void nfft_buf_t<Rank>::do_nfft() {
    if (nfft_type == nfft_type_t::type1)
      do_nfft_type1();
    else if (nfft_type == nfft_type_t::type3)
      do_nfft_type3();
    else
      do_direct();
  }

  // The static_assert on Rank in nfft_buf_t makes this exhaustive.
  template struct nfft_buf_t<1>;
  template struct nfft_buf_t<2>;
  template struct nfft_buf_t<3>;

} // namespace triqs::utility
