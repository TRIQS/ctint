#pragma once
#include <xsimd/xsimd.hpp>
#include <poet/poet.hpp>

#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <utility>

#include "../types.hpp"

// C99 restrict is not standard C++; every compiler triqs builds with spells it differently.
#if defined(_MSC_VER)
#define TRIQS_RESTRICT __restrict
#elif defined(__GNUC__) // also clang, icx, nvc++
#define TRIQS_RESTRICT __restrict__
#else
#define TRIQS_RESTRICT
#endif

// M3/M4 iw-accumulate kernels.
//
// Every measure reduces to the same block update, applied once per frequency-mesh point:
//
//   acc(i,j,k,l) += sign * M1a(j,i) * M2a(l,k)          (all block pairs)
//                 - sign * M1b(l,i) * M2b(j,k)          (diagonal block pairs only)
//
// Memory contract, enforced by static_assert in accumulate_block:
//   * the M-like operands store their target matrix transposed (nda stride order
//     {..., rank-1, rank-2}), so M2a(l,k) is M2a.data()[k*n + l] and M1b(:,i) is contiguous.
//     That is what makes l the contiguous axis of both terms, so each one scales a contiguous
//     run of complex by a coefficient that stays constant over the run.
//   * acc is C-contiguous, so its (i,j) planes are consecutive n*n runs.
//
// poet::dispatch turns the runtime block size into a compile-time N, so the accumulation loops
// have statically known trip counts and fully unroll. N == 0 selects the runtime-length path.
namespace triqs_ctint::measures {
  namespace {
    using cplx = std::complex<double>;

    // Largest block size served by a compile-time-length kernel.
    inline constexpr int max_block = 8;

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
    // assumed: an ISA need not provide a 128-bit double batch, so min_lanes is not always 2.
    inline constexpr long min_lanes = long(narrowest_batch<double, 1>());
    inline constexpr long max_lanes = long(widest_batch<double, 64>());

    template <long W> using vec = xsimd::make_sized_batch_t<double, W>;

    // --------------------------------------------- complex arithmetic on packed pairs

    // Complex data is interleaved (re,im,re,im,...), so a W-lane batch holds W/2 complex.
    // Swap real/imag within each pair: lane i <-> i^1.
    struct swap_re_im {
      static constexpr std::uint64_t get(std::uint64_t i, std::uint64_t) { return i ^ 1u; }
    };
    template <long W>
    inline constexpr auto swap_re_im_mask = xsimd::make_batch_constant<std::uint64_t, swap_re_im, typename vec<W>::arch_type>();

    // (cr + i*ci) * x: one swizzle + one vfmaddsub231pd.
    template <long W> XSIMD_INLINE vec<W> scale_by_complex(const vec<W> cr, const vec<W> ci, const vec<W> x) noexcept {
      return xsimd::fmas(cr, x, ci * xsimd::swizzle(x, swap_re_im_mask<W>));
    }

    // ------------------------------------------------------------------ length cover

    // Widest-first cover of D contiguous doubles: as many W-wide ops as fit, then halve W.
    // Halving stops at min_lanes, so a remainder narrower than that would be dropped silently.
    //
    // Each tile writes exactly the bytes it covers. Covering D with masked max_lanes tiles instead
    // issues 2.3x fewer instructions but runs 2.2x slower at N=3: a masked store writes a subset of
    // a max_lanes range that the next row's load overlaps, blocking store-to-load forwarding.
    template <long D, long W = max_lanes, long Off = 0, class Op> XSIMD_INLINE void for_each_tile(Op op) noexcept {
      static_assert(D % min_lanes == 0, "for_each_tile would leave an uncovered tail: D must be a multiple of min_lanes");
      if constexpr (Off >= D) {
        return;
      } else if constexpr (D - Off >= W) {
        op.template operator()<Off, W>();
        for_each_tile<D, W, Off + W>(op);
      } else if constexpr (W > min_lanes) {
        for_each_tile<D, W / 2, Off>(op);
      }
    }

    // -------------------------------------------------- scaled accumulation primitives
    //
    // Two operations: acc += c*x, and the fused acc += c1*x1 - c2*x2 in one pass over acc.
    // The suffix names the operand shape, not the arithmetic.

    // An (i,j)-invariant operand held in vector registers: 2*Len doubles as full_count widest
    // tiles plus at most one min_lanes tail. Len == N*N, so 2*Len % max_lanes is 0 (N even) or
    // min_lanes (N odd). Tile indices must be compile-time (poet::static_for) or the array
    // spills to memory. Len == 0 is the "nothing to hold" case and costs nothing.
    //
    // At the top of the range (Len == 64, N == 8) this wants 16 of the 32 zmm registers live
    // across the whole (i,j) loop and so spills: it trades spill stores for not re-loading the
    // operand on every (i,j).
    template <long Len> struct held_operand {
      static constexpr long lanes = 2 * Len, full_count = lanes / max_lanes, tail_lanes = lanes % max_lanes;
      static_assert(tail_lanes == 0 || tail_lanes == min_lanes);
      std::array<vec<max_lanes>, full_count> full;
      vec<min_lanes> tail;
      explicit held_operand(const cplx *TRIQS_RESTRICT x) noexcept {
        const auto *TRIQS_RESTRICT xd = reinterpret_cast<const double *>(x);
        poet::static_for<full_count>([&]<auto Tile>() { full[Tile] = vec<max_lanes>::load_unaligned(xd + Tile * max_lanes); });
        if constexpr (tail_lanes) tail = vec<min_lanes>::load_unaligned(xd + full_count * max_lanes);
      }
    };

    // acc[0:Len] += c * x, with x already held in registers.
    template <long Len> XSIMD_INLINE void add_scaled_held(const cplx c, const held_operand<Len> &x, cplx *TRIQS_RESTRICT acc) noexcept {
      using held              = held_operand<Len>;
      auto *TRIQS_RESTRICT ad = reinterpret_cast<double *>(acc);
      poet::static_for<held::full_count>([&]<auto Tile>() {
        using batch = vec<max_lanes>;
        auto av     = batch::load_unaligned(ad + Tile * max_lanes);
        av += scale_by_complex<max_lanes>(batch(c.real()), batch(c.imag()), x.full[Tile]);
        av.store_unaligned(ad + Tile * max_lanes);
      });
      if constexpr (held::tail_lanes) {
        using batch        = vec<min_lanes>;
        constexpr long off = held::full_count * max_lanes;
        auto av            = batch::load_unaligned(ad + off);
        av += scale_by_complex<min_lanes>(batch(c.real()), batch(c.imag()), x.tail);
        av.store_unaligned(ad + off);
      }
    }

    // acc[0:Len] += c1 * x1[0:Len] - c2 * x2[0:Len], both terms in one pass over acc.
    template <long Len>
    void add_scaled_minus_tiled(const cplx c1, const cplx *TRIQS_RESTRICT x1, const cplx c2, const cplx *TRIQS_RESTRICT x2,
                        cplx *TRIQS_RESTRICT acc) noexcept {
      const auto *TRIQS_RESTRICT x1d = reinterpret_cast<const double *>(x1);
      const auto *TRIQS_RESTRICT x2d = reinterpret_cast<const double *>(x2);
      auto *TRIQS_RESTRICT ad        = reinterpret_cast<double *>(acc);
      for_each_tile<2 * Len>([&]<long Off, long W>() {
        using batch = vec<W>;
        auto av     = batch::load_unaligned(ad + Off);
        av += scale_by_complex<W>(batch(c1.real()), batch(c1.imag()), batch::load_unaligned(x1d + Off))
           - scale_by_complex<W>(batch(c2.real()), batch(c2.imag()), batch::load_unaligned(x2d + Off));
        av.store_unaligned(ad + Off);
      });
    }

    // Small-block bodies for !prefer_simd. These stay separate functions: written as loops at
    // their call sites in accumulate_block they cost 12% at N=2.
    template <long Len> void add_scaled_plain(const cplx c, const cplx *TRIQS_RESTRICT x, cplx *TRIQS_RESTRICT acc) noexcept {
      for (long l = 0; l < Len; ++l) acc[l] += c * x[l];
    }
    template <long Len>
    void add_scaled_minus_plain(const cplx c1, const cplx *TRIQS_RESTRICT x1, const cplx c2, const cplx *TRIQS_RESTRICT x2,
                                cplx *TRIQS_RESTRICT acc) noexcept {
      for (long l = 0; l < Len; ++l) acc[l] += c1 * x1[l] - c2 * x2[l];
    }

    // Runtime-length variants, for blocks wider than max_block.
    void add_scaled_runtime(const cplx c, const cplx *TRIQS_RESTRICT x, cplx *TRIQS_RESTRICT acc, const long len) noexcept {
      const auto *TRIQS_RESTRICT xd = reinterpret_cast<const double *>(x);
      auto *TRIQS_RESTRICT ad       = reinterpret_cast<double *>(acc);
      using batch                   = vec<max_lanes>;
      const batch cr(c.real()), ci(c.imag());
      // max_lanes is a power of two, so the mask truncates to the last whole batch.
      const long vec_lanes = (2 * len) & -max_lanes;
      for (long off = 0; off < vec_lanes; off += max_lanes) {
        auto av = batch::load_unaligned(ad + off);
        av += scale_by_complex<max_lanes>(cr, ci, batch::load_unaligned(xd + off));
        av.store_unaligned(ad + off);
      }
      for (long l = vec_lanes / 2; l < len; ++l) acc[l] += c * x[l];
    }
    void add_scaled_minus_runtime(const cplx c1, const cplx *TRIQS_RESTRICT x1, const cplx c2, const cplx *TRIQS_RESTRICT x2, cplx *TRIQS_RESTRICT acc,
                          const long len) noexcept {
      const auto *TRIQS_RESTRICT x1d = reinterpret_cast<const double *>(x1);
      const auto *TRIQS_RESTRICT x2d = reinterpret_cast<const double *>(x2);
      auto *TRIQS_RESTRICT ad        = reinterpret_cast<double *>(acc);
      using batch                    = vec<max_lanes>;
      const batch c1r(c1.real()), c1i(c1.imag()), c2r(c2.real()), c2i(c2.imag());
      const long vec_lanes = (2 * len) & -max_lanes;
      for (long off = 0; off < vec_lanes; off += max_lanes) {
        auto av = batch::load_unaligned(ad + off);
        av += scale_by_complex<max_lanes>(c1r, c1i, batch::load_unaligned(x1d + off))
           - scale_by_complex<max_lanes>(c2r, c2i, batch::load_unaligned(x2d + off));
        av.store_unaligned(ad + off);
      }
      for (long l = vec_lanes / 2; l < len; ++l) acc[l] += c1 * x1[l] - c2 * x2[l];
    }

    // ------------------------------------------------------------ single-orbital rank-1 update

    // At bl_size == 1 every block axis has length 1, so accumulate_block has nothing to vectorize
    // and would pay its dispatch once per mesh point. Both full-mesh M3 updates are separable in
    // the two mesh indices there, so they accumulate along the contiguous iw2 axis instead: one
    // scaled add per acc row. Only the two full-mesh kernels are separable -- dlr2d's mesh is
    // sparse (no product structure) and the uniform iw3pp reads GM2 skewed as GM2[iW - iw].
    //
    // inline, not static: only the full-mesh kernels instantiate these, so internal linkage
    // would trip -Wunused-function in every other TU that includes this header.
    inline void add_scaled_outer(const cplx c, const cplx *TRIQS_RESTRICT x1, const cplx *TRIQS_RESTRICT x2, cplx *TRIQS_RESTRICT acc,
                            const long n1, const long n2) noexcept {
      for (long i1 = 0; i1 < n1; ++i1) add_scaled_runtime(c * x1[i1], x2, acc + i1 * n2, n2);
    }

    // acc(iw1, iw2) += c * m(iw1, iw2) - c2 * x1(iw1) * x2(iw2), both terms in one pass over acc.
    inline void add_scaled_minus_outer(const cplx c, const cplx *TRIQS_RESTRICT m, const cplx c2, const cplx *TRIQS_RESTRICT x1,
                                     const cplx *TRIQS_RESTRICT x2, cplx *TRIQS_RESTRICT acc, const long n1, const long n2) noexcept {
      for (long i1 = 0; i1 < n1; ++i1) add_scaled_minus_runtime(c, m + i1 * n2, c2 * x1[i1], x2, acc + i1 * n2, n2);
    }

    // ------------------------------------------------------------ block accumulation

    // Hand-vectorize over (k,l) only for blocks spanning at least two widest vectors: N >= 3 on
    // AVX-512, N >= 2 on AVX2/SSE. Below that the plain loop is faster, even though clang and gcc
    // auto-vectorize it to the same swizzle+vfmaddsub form -- the cost sits in the surrounding
    // view machinery rather than the vector ops, so the threshold is empirical.
    //
    // It keys on the block size N, not on the individual run length: the diagonal path gates an
    // N-long fused add on it, so at N=3 on AVX-512 a 3-complex run takes the vector body even
    // though it is shorter than one widest vector.
    template <int N> inline constexpr bool prefer_simd = long{N} * N >= max_lanes;

    // acc(i,j,k,l) += sign * M1a(j,i) * M2a(l,k)  [- sign * M1b(l,i) * M2b(j,k) on the diagonal].
    // N == bl2_size at compile time; N == 0 selects the runtime-length path.
    template <int N, bool diagonal>
    void accumulate_block(const mc_weight_t sign, const auto &M1a, const auto &M2a, const auto &M1b, const auto &M2b, auto &acc,
                          const long bl1_size, const long bl2_size) noexcept {
      static_assert(std::decay_t<decltype(M2a.indexmap())>::is_stride_order_Fortran(),
                    "M2a must store its target transposed: the kernel reads M2a(l,k) as M2a.data()[k*n + l]");
      static_assert(std::decay_t<decltype(M1b.indexmap())>::is_stride_order_Fortran(),
                    "M1b must store its target transposed: the kernel reads the column M1b(:,i) contiguously");
      static_assert(std::decay_t<decltype(acc.indexmap())>::is_stride_order_C(), "acc must be C-ordered: its (i,j) planes are walked linearly");
      assert(acc.indexmap().is_contiguous());

      // M2a is (i,j)-invariant: hold it in registers once, not once per (i,j). Length 0 (every
      // path that does not want it) makes held_operand empty and its construction a no-op.
      constexpr long held_len                   = (N > 0 && !diagonal && prefer_simd<N>) ? long{N} * N : 0;
      const cplx *const TRIQS_RESTRICT M2a_flat = M2a.data();
      const held_operand<held_len> M2a_held{M2a_flat};

      // The acc(i,j,:,:) planes are contiguous and consecutive, so walk them with a pointer.
      // Re-deriving &acc(i,j,0,0) costs an nda 4-index offset plus view refcounting, which at
      // small N dominates the arithmetic.
      const long plane_len       = (N > 0) ? long{N} * N : bl2_size * bl2_size;
      cplx *TRIQS_RESTRICT plane = acc.data();
      for (long i = 0; i < bl1_size; ++i) {
        const cplx *const TRIQS_RESTRICT M1b_col = &M1b(0, i);
        for (long j = 0; j < bl1_size; ++j, plane += plane_len) {
          const cplx c1 = M1a(j, i) * sign;
          if constexpr (N == 0) {
            if constexpr (!diagonal) {
              add_scaled_runtime(c1, M2a_flat, plane, plane_len);
            } else {
              for (long k = 0; k < bl2_size; ++k)
                add_scaled_minus_runtime(c1, M2a_flat + k * bl2_size, M2b(j, k) * sign, M1b_col, plane + k * bl2_size, bl2_size);
            }
          } else if constexpr (!diagonal) {
            if constexpr (prefer_simd<N>)
              add_scaled_held<held_len>(c1, M2a_held, plane);
            else
              add_scaled_plain<long{N} * N>(c1, M2a_flat, plane);
          } else {
            for (int k = 0; k < N; ++k) {
              const cplx c2 = M2b(j, k) * sign;
              if constexpr (prefer_simd<N>)
                add_scaled_minus_tiled<N>(c1, M2a_flat + long{k} * N, c2, M1b_col, plane + long{k} * N);
              else
                add_scaled_minus_plain<N>(c1, M2a_flat + long{k} * N, c2, M1b_col, plane + long{k} * N);
            }
          }
        }
      }
    }

    // The three M3 pp kernels below pass the same operand pair twice, as (M1a, M2a, M1a, M2a).
    // At N == 1 that makes both terms of the diagonal update M1a(0,0) * M2a(0,0), so they cancel
    // exactly and the whole block pair contributes nothing -- half the block pairs of a
    // single-orbital model. N > 1 does not cancel (the two terms read different index pairs), and
    // neither does any ph kernel or iw4pp (distinct M1b/M2b).
    template <int N, bool diagonal> inline constexpr bool pp_pair_cancels = (N == 1 && diagonal);

    // ======================== M4 kernels (uniform mesh) ========================

    constexpr auto iw4_accumulate_kernel = []<int N, bool diagonal>(const mc_weight_t sign, const auto &M, auto &M4_iw,
                                                                    const auto bl1, const auto bl2) noexcept {
      auto &[iw_mesh, _, _] = M4_iw(0, 0).mesh();
      auto const &M1        = M[bl1];
      auto const &M2        = M[bl2];
      auto const bl1_size   = M1.target_shape()[0];
      auto const bl2_size   = M2.target_shape()[0];
      auto &M4              = M4_iw(bl1, bl2);

      for (const auto &iw1 : iw_mesh) {
        for (const auto &iw2 : iw_mesh) {
          for (const auto &iw3 : iw_mesh) {
            const auto iw4 = iw1 + iw3 - iw2;
            const auto M1a = M1[iw2.value(), iw1];
            const auto M2a = M2[iw4, iw3];
            const auto M1b = M1[iw4, iw1];
            const auto M2b = M2[iw2.value(), iw3];
            auto acc       = M4[iw1, iw2, iw3];
            accumulate_block<N, diagonal>(sign, M1a, M2a, M1b, M2b, acc, bl1_size, bl2_size);
          }
        }
      }
    };

    constexpr auto iw4ph_accumulate_kernel = []<int N, bool diagonal>(const mc_weight_t sign, const auto &M, auto &M4_iw,
                                                                      const auto bl1, const auto bl2) noexcept {
      auto const &[iW_mesh, iw_mesh, _] = M4_iw(0, 0).mesh();
      auto const &M1                    = M[bl1];
      auto const &M2                    = M[bl2];
      auto const bl1_size               = M1.target_shape()[0];
      auto const bl2_size               = M2.target_shape()[0];
      auto &M4                          = M4_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          for (auto iwp : iw_mesh) {
            const auto M1a = M1[iW + iw, iw.value()];
            const auto M2a = M2[iwp.value(), iW + iwp];
            const auto M1b = M1[iwp.value(), iw.value()];
            const auto M2b = M2[iW + iw, iW + iwp];
            auto acc       = M4[iW, iw, iwp];
            accumulate_block<N, diagonal>(sign, M1a, M2a, M1b, M2b, acc, bl1_size, bl2_size);
          }
        }
      }
    };

    constexpr auto iw4pp_accumulate_kernel = []<int N, bool diagonal>(const mc_weight_t sign, const auto &M, auto &M4_iw,
                                                                      const auto bl1, const auto bl2) noexcept {
      auto const &[iW_mesh, iw_mesh, _] = M4_iw(0, 0).mesh();
      auto const &M1                    = M[bl1];
      auto const &M2                    = M[bl2];
      auto const bl1_size               = M1.target_shape()[0];
      auto const bl2_size               = M2.target_shape()[0];
      auto &M4                          = M4_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          for (auto iwp : iw_mesh) {
            const auto M1a = M1[iW - iwp, iw.value()];
            const auto M2a = M2[iwp.value(), iW - iw];
            const auto M1b = M1[iwp.value(), iw.value()];
            const auto M2b = M2[iW - iwp, iW - iw];
            auto acc       = M4[iW, iw, iwp];
            accumulate_block<N, diagonal>(sign, M1a, M2a, M1b, M2b, acc, bl1_size, bl2_size);
          }
        }
      }
    };

    // ======================== M3 kernels (uniform mesh) ========================

    constexpr auto iw3pp_accumulate_kernel = []<int N, bool diagonal>(const mc_weight_t sign, const auto &GM, auto &M3_iw,
                                                                      const auto bl1, const auto bl2) noexcept {
      if constexpr (pp_pair_cancels<N, diagonal>) return;
      auto const [iW_mesh, iw_mesh] = M3_iw(0, 0).mesh();
      auto const bl1_size           = GM[bl1].target_shape()[0];
      auto const bl2_size           = GM[bl2].target_shape()[0];
      auto const &GM1               = GM[bl1];
      auto const &GM2               = GM[bl2];
      auto &M3                      = M3_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          const auto M1a = GM1[iw.value()];
          const auto M2a = GM2[iW - iw];
          auto acc       = M3[iW, iw];
          accumulate_block<N, diagonal>(sign, M1a, M2a, M1a, M2a, acc, bl1_size, bl2_size);
        }
      }
    };

    constexpr auto iw3ph_accumulate_kernel = []<int N, bool diagonal>(const mc_weight_t sign, const auto &M, const auto &GMG,
                                                                      const auto &GM, const auto &MG, auto &M3_iw,
                                                                      const auto bl1, const auto bl2) noexcept {
      auto const [iW_mesh, iw_mesh] = M3_iw(0, 0).mesh();
      auto const bl1_size           = M[bl1].target_shape()[0];
      auto const bl2_size           = M[bl2].target_shape()[0];
      auto const &M1                = M[bl1];
      auto const &GMG2              = GMG(bl2);
      auto const &GM1               = GM[bl1];
      auto const &MG2               = MG(bl2);
      auto &M3                      = M3_iw(bl1, bl2);

      for (auto iW : iW_mesh) {
        for (auto iw : iw_mesh) {
          const auto M1a = M1[iW + iw, iw.value()];
          const auto M1b = GM1[iw.value()];
          const auto M2b = MG2[iW + iw];
          auto acc       = M3[iW, iw];
          accumulate_block<N, diagonal>(sign, M1a, GMG2, M1b, M2b, acc, bl1_size, bl2_size);
        }
      }
    };

    // ======================== M3 DLR2D kernels ========================

    constexpr auto dlr2d_iw3ph_accumulate_kernel =
       []<int N, bool diagonal>(const mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG,
                                auto &M3_iw, const auto bl1, const auto bl2) noexcept {
         auto const bl1_size = M[bl1].target_shape()[0];
         auto const bl2_size = GMG(bl2).shape()[0];
         auto const &M1      = M[bl1];
         auto const &GMG2    = GMG(bl2);
         auto const &GM1     = GM[bl1];
         auto const &MG2     = MG[bl2];
         auto &M3            = M3_iw(bl1, bl2);

         for (auto mp : M3.mesh()) {
           auto [iw1, iw2] = mp.value();
           const auto M1a  = M1[mp];
           const auto M1b  = GM1[iw1];
           const auto M2b  = MG2[iw2];
           auto acc        = M3[mp];
           accumulate_block<N, diagonal>(sign, M1a, GMG2, M1b, M2b, acc, bl1_size, bl2_size);
         }
       };

    constexpr auto dlr2d_iw3pp_accumulate_kernel = []<int N, bool diagonal>(const mc_weight_t sign, const auto &GM,
                                                                            auto &M3_iw, const auto bl1,
                                                                            const auto bl2) noexcept {
      if constexpr (pp_pair_cancels<N, diagonal>) return;
      auto const bl1_size = GM[bl1].target_shape()[0];
      auto const bl2_size = GM[bl2].target_shape()[0];
      auto const &GM1     = GM[bl1];
      auto const &GM2     = GM[bl2];
      auto &M3            = M3_iw(bl1, bl2);

      for (auto mp : M3.mesh()) {
        auto [iw1, iw2] = mp.value();
        const auto M1a  = GM1[iw1];
        const auto M2a  = GM2[iw2];
        auto acc        = M3[mp];
        accumulate_block<N, diagonal>(sign, M1a, M2a, M1a, M2a, acc, bl1_size, bl2_size);
      }
    };

    // ======================== M3 full (uniform prod mesh) kernels ========================

    constexpr auto full_iw3ph_accumulate_kernel =
       []<int N, bool diagonal>(const mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG,
                                auto &M3_iw, const auto bl1, const auto bl2) noexcept {
         auto const bl1_size = M[bl1].target_shape()[0];
         auto const bl2_size = GMG(bl2).shape()[0];
         auto const &M1      = M[bl1];
         auto const &GMG2    = GMG(bl2);
         auto const &GM1     = GM[bl1];
         auto const &MG2     = MG[bl2];
         auto &M3            = M3_iw(bl1, bl2);

         // acc(iw1, iw2) += sign * GMG2 * M(iw1, iw2) [- sign * GM1(iw1) * MG2(iw2) on the diagonal].
         // Term 1 has a constant coefficient over the whole plane and M shares M3's mesh, so off the
         // diagonal it is one scaled add over the full plane; on the diagonal the two terms fuse per row.
         if constexpr (N == 1) {
           const auto n1 = M3.data().shape()[0], n2 = M3.data().shape()[1];
           const auto c1 = sign * GMG2(0, 0);
           if constexpr (diagonal)
             add_scaled_minus_outer(c1, M1.data().data(), sign, GM1.data().data(), MG2.data().data(), M3.data().data(), n1, n2);
           else
             add_scaled_runtime(c1, M1.data().data(), M3.data().data(), n1 * n2);
           return;
         }

         for (auto mp : M3.mesh()) {
           auto [mp1, mp2] = mp;
           auto iw1        = mp1.value();
           auto iw2        = mp2.value();
           const auto M1a  = M1[closest_mesh_pt(iw1, iw2)];
           const auto M1b  = GM1[iw1];
           const auto M2b  = MG2[iw2];
           auto acc        = M3[mp];
           accumulate_block<N, diagonal>(sign, M1a, GMG2, M1b, M2b, acc, bl1_size, bl2_size);
         }
       };

    constexpr auto full_iw3pp_accumulate_kernel = []<int N, bool diagonal>(const mc_weight_t sign, const auto &GM,
                                                                           auto &M3_iw, const auto bl1,
                                                                           const auto bl2) noexcept {
      if constexpr (pp_pair_cancels<N, diagonal>) return;
      auto const bl1_size = GM[bl1].target_shape()[0];
      auto const bl2_size = GM[bl2].target_shape()[0];
      auto const &GM1     = GM[bl1];
      auto const &GM2     = GM[bl2];
      auto &M3            = M3_iw(bl1, bl2);

      // acc(iw1, iw2) += sign * GM1(iw1) * GM2(iw2). Diagonal pairs never reach here (they cancel).
      if constexpr (N == 1) {
        add_scaled_outer(sign, GM1.data().data(), GM2.data().data(), M3.data().data(), M3.data().shape()[0], M3.data().shape()[1]);
        return;
      }

      for (auto mp : M3.mesh()) {
        auto [mp1, mp2] = mp;
        auto iw1        = mp1.value();
        auto iw2        = mp2.value();
        const auto M1a  = GM1[iw1];
        const auto M2a  = GM2[iw2];
        auto acc        = M3[mp];
        accumulate_block<N, diagonal>(sign, M1a, M2a, M1a, M2a, acc, bl1_size, bl2_size);
      }
    };
  } // namespace

  namespace simd {

    // Route the runtime (bl2_size, diagonal) pair to a compile-time kernel specialization.
    // bl2_size outside [1, max_block] maps to N == 0, the runtime-length path.
    template <auto kernel> void kernel_dispatch(const auto bl2_size, const bool diagonal, auto &&...args) noexcept {
      const int n = (bl2_size <= max_block) ? static_cast<int>(bl2_size) : 0;
      poet::dispatch([&]<int N, int D>() { kernel.template operator()<N, static_cast<bool>(D)>(std::forward<decltype(args)>(args)...); },
                     poet::dispatch_param<poet::inclusive_range<0, max_block>>{n},
                     poet::dispatch_param<poet::inclusive_range<0, 1>>{diagonal});
    }

    // M4 dispatchers
    void iw4_accumulate(const mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) noexcept {
      kernel_dispatch<iw4_accumulate_kernel>(bl2_size, bl1 == bl2, sign, M, M4_iw, bl1, bl2);
    }

    void iw4ph_accumulate(const mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) noexcept {
      kernel_dispatch<iw4ph_accumulate_kernel>(bl2_size, bl1 == bl2, sign, M, M4_iw, bl1, bl2);
    }

    void iw4pp_accumulate(const mc_weight_t sign, const auto &M, auto &M4_iw, const auto bl1, const auto bl2, const auto bl2_size) noexcept {
      kernel_dispatch<iw4pp_accumulate_kernel>(bl2_size, bl1 == bl2, sign, M, M4_iw, bl1, bl2);
    }

    // M3 uniform mesh dispatchers
    void iw3ph_accumulate(const mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG, auto &M3_iw, const auto bl1,
                          const auto bl2, const auto bl2_size) noexcept {
      kernel_dispatch<iw3ph_accumulate_kernel>(bl2_size, bl1 == bl2, sign, M, GMG, GM, MG, M3_iw, bl1, bl2);
    }

    void iw3pp_accumulate(const mc_weight_t sign, const auto &GM, auto &M3_iw, const auto bl1, const auto bl2, const auto bl2_size) noexcept {
      kernel_dispatch<iw3pp_accumulate_kernel>(bl2_size, bl1 == bl2, sign, GM, M3_iw, bl1, bl2);
    }

    // M3 DLR2D dispatchers
    void dlr2d_iw3ph_accumulate(const mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG, auto &M3_iw,
                                const auto bl1, const auto bl2, const auto bl2_size) noexcept {
      kernel_dispatch<dlr2d_iw3ph_accumulate_kernel>(bl2_size, bl1 == bl2, sign, M, GMG, GM, MG, M3_iw, bl1, bl2);
    }

    void dlr2d_iw3pp_accumulate(const mc_weight_t sign, const auto &GM, auto &M3_iw, const auto bl1, const auto bl2,
                                const auto bl2_size) noexcept {
      kernel_dispatch<dlr2d_iw3pp_accumulate_kernel>(bl2_size, bl1 == bl2, sign, GM, M3_iw, bl1, bl2);
    }

    // M3 full (uniform prod mesh) dispatchers
    void full_iw3ph_accumulate(const mc_weight_t sign, const auto &M, const auto &GMG, const auto &GM, const auto &MG, auto &M3_iw,
                               const auto bl1, const auto bl2, const auto bl2_size) noexcept {
      kernel_dispatch<full_iw3ph_accumulate_kernel>(bl2_size, bl1 == bl2, sign, M, GMG, GM, MG, M3_iw, bl1, bl2);
    }

    void full_iw3pp_accumulate(const mc_weight_t sign, const auto &GM, auto &M3_iw, const auto bl1, const auto bl2,
                               const auto bl2_size) noexcept {
      kernel_dispatch<full_iw3pp_accumulate_kernel>(bl2_size, bl1 == bl2, sign, GM, M3_iw, bl1, bl2);
    }

  } // namespace simd
} // namespace triqs_ctint::measures
