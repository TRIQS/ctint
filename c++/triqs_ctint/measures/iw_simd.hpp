// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <xsimd/xsimd.hpp>

#include <complex>
#include <cstdint>
#include <type_traits>

// SIMD machinery shared by the M3/M4 iw-accumulate loops. Only the primitives live here; each
// measure writes its own accumulation loop over them.

// C99 restrict is not standard C++; every compiler triqs builds with spells it differently.
#if defined(_MSC_VER)
#define TRIQS_RESTRICT __restrict
#elif defined(__GNUC__) // also clang, icx, nvc++
#define TRIQS_RESTRICT __restrict__
#else
#define TRIQS_RESTRICT
#endif

namespace triqs_ctint::measures {
  namespace {

    // ------------------------------------------------------------------- lane widths

    // Narrowest/widest double batch the ISA actually provides.
    template <class T, auto N> constexpr auto narrowest_batch() {
      if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, N>>) {
        return narrowest_batch<T, N * 2>();
      } else {
        return N;
      }
    }
    template <class T, auto N> constexpr auto widest_batch() {
      if constexpr (std::is_void_v<xsimd::make_sized_batch_t<T, N>>) {
        return widest_batch<T, N / 2>();
      } else {
        return N;
      }
    }

    // Counted in doubles: 2 -> 128-bit, 4 -> 256-bit, 8 -> 512-bit. Both are discovered, not
    // assumed: an ISA need not provide a 128-bit double batch, so min_simd is not always 2.
    inline constexpr long min_simd = long(narrowest_batch<double, 1>());
    inline constexpr long max_simd = long(widest_batch<double, 64>());

    template <long W> using vec = xsimd::make_sized_batch_t<double, W>;

    // --------------------------------------------- complex arithmetic on packed pairs

    // Complex data is interleaved (re,im,re,im,...), so a W-lane batch holds W/2 complex.
    // Swap real/imag within each pair: lane i <-> i^1.
    struct swap_re_im {
      static constexpr std::uint64_t get(std::uint64_t i, std::uint64_t) { return i ^ 1u; }
    };
    template <class B> inline constexpr auto swap_re_im_mask = xsimd::make_batch_constant<std::uint64_t, swap_re_im, typename B::arch_type>();

    // (cr + i*ci) * x, both packed: one swizzle + one vfmaddsub231pd.
    template <class B> XSIMD_INLINE B cmul(const B cr, const B ci, const B x) noexcept {
      return xsimd::fmas(cr, x, ci * xsimd::swizzle(x, swap_re_im_mask<B>));
    }

    // acc + (cr + i*ci) * x. The add cannot fold into the fmaddsub -- that one already spends its
    // addend on the alternating-sign term -- so this is cmul plus one vaddpd.
    template <class B> XSIMD_INLINE B cfma(const B cr, const B ci, const B x, const B acc) noexcept { return acc + cmul(cr, ci, x); }

    // ------------------------------------------------------------------------ block sizes

    using cplx = std::complex<double>;

    // Largest block size served by a compile-time-length kernel. Above it the accumulation loops
    // take their runtime-length path.
    inline constexpr int max_block = 8;

    // True once an N*N plane's 2*N*N doubles fill a widest vector, which is what makes it worth
    // keeping that operand live in registers across the (i,j) loops. Below it the operand is
    // narrower than one vector and fma_run reloads it, which costs nothing extra.
    //
    // This decides register residency only. Every run is vectorized whatever N is: fma_run
    // walks the power-of-two widths from max_simd down to min_simd and takes the widest that
    // fits, so N == 1 gets a 2-double batch rather than a scalar loop.
    template <int N> inline constexpr bool worth_holding = 2 * long{N} * N >= max_simd;

  } // namespace
} // namespace triqs_ctint::measures
