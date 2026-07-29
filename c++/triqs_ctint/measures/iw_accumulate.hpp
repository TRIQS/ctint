#pragma once
#include <poet/poet.hpp>

#include <array>
#include <cassert>
#include <complex>
#include <utility>

#include "../types.hpp"
#include "./iw_simd.hpp"

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
    // ------------------------------------------------------------ block accumulation

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

      const cplx *const TRIQS_RESTRICT M2a_flat = M2a.data();
      const auto *const TRIQS_RESTRICT M2a_d    = reinterpret_cast<const double *>(M2a_flat);

      // The !diagonal path reads all N*N of M2a on every (i,j), so hold it in vector registers
      // for the whole loop: 2*N*N doubles as `held_full` widest tiles plus, when N is odd, one
      // min_simd tail. Tile indices must be compile-time (poet::static_for) or the array spills.
      //
      // At N == 8 this wants 16 of the 32 zmm live across the loop and does spill: it trades
      // spill stores for not re-loading the operand on every (i,j).
      constexpr bool hold      = N > 0 && !diagonal && prefer_simd<N>;
      constexpr long held_full = hold ? 2 * long{N} * N / max_simd : 0;
      constexpr long held_tail = hold ? 2 * long{N} * N % max_simd : 0;
      static_assert(held_tail == 0 || held_tail == min_simd);
      std::array<vec<max_simd>, held_full> M2a_full;
      vec<min_simd> M2a_tail{};
      if constexpr (hold) {
        const auto *TRIQS_RESTRICT xd = reinterpret_cast<const double *>(M2a_flat);
        poet::static_for<held_full>([&]<auto Tile>() { M2a_full[Tile] = vec<max_simd>::load_unaligned(xd + Tile * max_simd); });
        if constexpr (held_tail) M2a_tail = vec<min_simd>::load_unaligned(xd + held_full * max_simd);
      }

      // The acc(i,j,:,:) planes are contiguous and consecutive, so walk them with a pointer.
      // Re-deriving &acc(i,j,0,0) costs an nda 4-index offset plus view refcounting, which at
      // small N dominates the arithmetic.
      const long plane_len       = (N > 0) ? long{N} * N : bl2_size * bl2_size;
      cplx *TRIQS_RESTRICT plane = acc.data();
      using wide                 = vec<max_simd>;
      for (long i = 0; i < bl1_size; ++i) {
        const cplx *const TRIQS_RESTRICT M1b_col = &M1b(0, i);
        // The restrict-qualified base pointers stay outside the k loops below. Declared inside one,
        // each iteration opens a fresh restrict scope, so clang can no longer tell that the store to
        // acc leaves M1b(:,i) alone and re-loads plus re-swizzles it on every k.
        const auto *const TRIQS_RESTRICT M1b_d = reinterpret_cast<const double *>(M1b_col);
        for (long j = 0; j < bl1_size; ++j, plane += plane_len) {
          const cplx c1                    = M1a(j, i) * sign;
          auto *const TRIQS_RESTRICT acc_d = reinterpret_cast<double *>(plane);

          if constexpr (N == 0) {
            // Runtime block size: whole widest vectors then a scalar tail. max_simd is a power of
            // two, so the mask truncates to the last whole batch.
            if constexpr (!diagonal) {
              const long lanes = (2 * plane_len) & -max_simd;
              const wide cr(c1.real()), ci(c1.imag());
              for (long off = 0; off < lanes; off += max_simd) {
                auto av = wide::load_unaligned(acc_d + off);
                av += cmul(cr, ci, wide::load_unaligned(M2a_d + off));
                av.store_unaligned(acc_d + off);
              }
              for (long l = lanes / 2; l < plane_len; ++l) plane[l] += c1 * M2a_flat[l];
            } else {
              const long lanes = (2 * bl2_size) & -max_simd;
              const wide c1r(c1.real()), c1i(c1.imag());
              for (long k = 0; k < bl2_size; ++k) {
                const cplx c2 = M2b(j, k) * sign;
                const wide c2r(c2.real()), c2i(c2.imag());
                const long base = 2 * k * bl2_size;
                for (long off = 0; off < lanes; off += max_simd) {
                  auto av = wide::load_unaligned(acc_d + base + off);
                  av += cmul(c1r, c1i, wide::load_unaligned(M2a_d + base + off)) - cmul(c2r, c2i, wide::load_unaligned(M1b_d + off));
                  av.store_unaligned(acc_d + base + off);
                }
                for (long l = lanes / 2; l < bl2_size; ++l) plane[k * bl2_size + l] += c1 * M2a_flat[k * bl2_size + l] - c2 * M1b_col[l];
              }
            }

          } else if constexpr (!diagonal) {
            // plane[0 : N*N] += c1 * M2a, the operand already in registers.
            if constexpr (prefer_simd<N>) {
              const wide cr(c1.real()), ci(c1.imag());
              poet::static_for<held_full>([&]<auto Tile>() {
                auto av = wide::load_unaligned(acc_d + Tile * max_simd);
                av += cmul(cr, ci, M2a_full[Tile]);
                av.store_unaligned(acc_d + Tile * max_simd);
              });
              if constexpr (held_tail) {
                using narrow = vec<min_simd>;
                auto av      = narrow::load_unaligned(acc_d + held_full * max_simd);
                av += cmul(narrow(c1.real()), narrow(c1.imag()), M2a_tail);
                av.store_unaligned(acc_d + held_full * max_simd);
              }
            } else
              for (long l = 0; l < long{N} * N; ++l) plane[l] += c1 * M2a_flat[l];

          } else {
            for (int k = 0; k < N; ++k) {
              const cplx c2 = M2b(j, k) * sign;
              // plane[k*N : k*N + N) += c1 * M2a(:,k) - c2 * M1b(:,i), both terms in one pass.
              if constexpr (prefer_simd<N>) {
                // 2*N contiguous doubles: as many widest ops as fit, then the even remainder as at
                // most one 4-wide plus one 2-wide op. Each op writes exactly the bytes it covers.
                // Covering the remainder with a masked widest store instead issues 2.3x fewer
                // instructions but runs 2.2x slower at N=3: that store writes a subset of a widest
                // range which the next k's load overlaps, blocking store-to-load forwarding.
                static_assert(min_simd == 2 && max_simd <= 8,
                              "the remainder cover below enumerates the cases for min_simd == 2, max_simd <= 8");
                constexpr long len = 2 * long{N}, rem = len % max_simd;
                const long base    = long{k} * len;
                auto op            = [&](auto w, const long off) {
                  using batch = typename decltype(w)::type;
                  auto av     = batch::load_unaligned(acc_d + base + off);
                  av += cmul(batch(c1.real()), batch(c1.imag()), batch::load_unaligned(M2a_d + base + off))
                     - cmul(batch(c2.real()), batch(c2.imag()), batch::load_unaligned(M1b_d + off));
                  av.store_unaligned(acc_d + base + off);
                };
                long off = 0;
                for (; off + max_simd <= len; off += max_simd) op(width<max_simd>, off);
                if constexpr (rem >= 4) {
                  op(width<4>, off);
                  off += 4;
                }
                if constexpr (rem % 4) op(width<min_simd>, off);
              } else
                for (int l = 0; l < N; ++l) plane[k * N + l] += c1 * M2a_flat[k * N + l] - c2 * M1b_col[l];
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
         // At bl_size == 1 every block axis has length 1, so accumulate_block has nothing to
         // vectorize and would pay its dispatch once per mesh point. Both terms are separable in
         // the two mesh indices here, so accumulate along the contiguous iw2 axis instead: one
         // scaled add per acc row. Only the two full-mesh kernels are separable -- dlr2d's mesh is
         // sparse (no product structure) and the uniform iw3pp reads GM2 skewed as GM2[iW - iw].
         if constexpr (N == 1) {
           const long n1                          = M3.data().shape()[0], n2 = M3.data().shape()[1];
           const cplx c1                          = sign * GMG2(0, 0);
           const cplx *const TRIQS_RESTRICT m     = M1.data().data();
           cplx *const TRIQS_RESTRICT acc         = M3.data().data();
           const auto *const TRIQS_RESTRICT m_d   = reinterpret_cast<const double *>(m);
           auto *const TRIQS_RESTRICT acc_d       = reinterpret_cast<double *>(acc);
           using wide                             = vec<max_simd>;
           const wide c1r(c1.real()), c1i(c1.imag());
           if constexpr (diagonal) {
             // acc(iw1,:) += c1 * M(iw1,:) - sign * GM1(iw1) * MG2(:), both terms in one pass.
             const cplx *const TRIQS_RESTRICT gm1   = GM1.data().data();
             const auto *const TRIQS_RESTRICT mg2_d = reinterpret_cast<const double *>(MG2.data().data());
             const cplx *const TRIQS_RESTRICT mg2   = MG2.data().data();
             const long lanes                       = (2 * n2) & -max_simd;
             for (long i1 = 0; i1 < n1; ++i1) {
               const cplx c2 = sign * gm1[i1];
               const wide c2r(c2.real()), c2i(c2.imag());
               const long base = 2 * i1 * n2;
               for (long off = 0; off < lanes; off += max_simd) {
                 auto av = wide::load_unaligned(acc_d + base + off);
                 av += cmul(c1r, c1i, wide::load_unaligned(m_d + base + off)) - cmul(c2r, c2i, wide::load_unaligned(mg2_d + off));
                 av.store_unaligned(acc_d + base + off);
               }
               for (long l = lanes / 2; l < n2; ++l) acc[i1 * n2 + l] += c1 * m[i1 * n2 + l] - c2 * mg2[l];
             }
           } else {
             // Term 1 has a constant coefficient over the whole plane and M shares M3's mesh, so it
             // is one scaled add over the flat plane.
             const long len = n1 * n2, lanes = (2 * len) & -max_simd;
             for (long off = 0; off < lanes; off += max_simd) {
               auto av = wide::load_unaligned(acc_d + off);
               av += cmul(c1r, c1i, wide::load_unaligned(m_d + off));
               av.store_unaligned(acc_d + off);
             }
             for (long l = lanes / 2; l < len; ++l) acc[l] += c1 * m[l];
           }
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
      // Same separable single-orbital shortcut as full_iw3ph: a rank-1 outer product, one scaled
      // add per acc row along the contiguous iw2 axis.
      if constexpr (N == 1) {
        const long n1                          = M3.data().shape()[0], n2 = M3.data().shape()[1];
        const cplx *const TRIQS_RESTRICT gm1   = GM1.data().data();
        const cplx *const TRIQS_RESTRICT gm2   = GM2.data().data();
        cplx *const TRIQS_RESTRICT acc         = M3.data().data();
        const auto *const TRIQS_RESTRICT gm2_d = reinterpret_cast<const double *>(gm2);
        auto *const TRIQS_RESTRICT acc_d       = reinterpret_cast<double *>(acc);
        using wide                             = vec<max_simd>;
        const long lanes                       = (2 * n2) & -max_simd;
        for (long i1 = 0; i1 < n1; ++i1) {
          const cplx c = sign * gm1[i1];
          const wide cr(c.real()), ci(c.imag());
          const long base = 2 * i1 * n2;
          for (long off = 0; off < lanes; off += max_simd) {
            auto av = wide::load_unaligned(acc_d + base + off);
            av += cmul(cr, ci, wide::load_unaligned(gm2_d + off));
            av.store_unaligned(acc_d + base + off);
          }
          for (long l = lanes / 2; l < n2; ++l) acc[i1 * n2 + l] += c * gm2[l];
        }
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
