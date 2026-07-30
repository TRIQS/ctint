// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3ph_iw_full.hpp"
#include "./iw_simd.hpp"

#include <poet/poet.hpp>

#include <array>
#include <cassert>

namespace triqs_ctint::measures {

  namespace {
    // One contiguous run of interleaved complex, len doubles == len/2 values:
    //
    //   acc[m] += c1 * a[m]                 (Diag == false)
    //   acc[m] += c1 * a[m] - c2 * b[m]     (Diag == true)
    //
    // Len == 0 defers the length to the runtime argument. Widest batch first, then half, down to
    // min_simd; whatever is narrower than the narrowest batch finishes one complex at a time. With
    // Len > 0 every bound here is a constant, so the loops unroll and nothing branches.
    //
    // Only acc is restrict: that is the one fact the loop needs (nothing else touches the
    // accumulator). a and b are read-only, and the non-diagonal calls pass the same pointer twice.
    template <bool Diag, long Len, long Width = max_simd>
    XSIMD_INLINE void fma_run(double *TRIQS_RESTRICT acc, const double *a, const double *b, const cplx c1, const cplx c2,
                              const long runtime_len = 0) noexcept {
      using batch        = vec<Width>;
      const long len     = (Len > 0) ? Len : runtime_len;
      const long covered = len - len % Width; // doubles this width takes; the rest drops to half
      const batch c1r(c1.real()), c1i(c1.imag());
      for (long start = 0; start < covered; start += Width) {
        auto av = cfma(c1r, c1i, batch::load_unaligned(a + start), batch::load_unaligned(acc + start));
        if constexpr (Diag) av -= cmul(batch(c2.real()), batch(c2.imag()), batch::load_unaligned(b + start));
        av.store_unaligned(acc + start);
      }

      constexpr long tail_len = Len % Width;
      constexpr long half     = Width / 2;
      if constexpr (Len == 0 || tail_len != 0) {
        if constexpr (half >= min_simd) {
          fma_run<Diag, tail_len, half>(acc + covered, a + covered, b + covered, c1, c2, len - covered);
        } else {
          const long rest                  = (len - covered) / 2; // leftover complex, narrower than any batch
          auto *const TRIQS_RESTRICT acc_c = reinterpret_cast<cplx *>(acc + covered);
          const auto *const a_c            = reinterpret_cast<const cplx *>(a + covered);
          const auto *const b_c            = reinterpret_cast<const cplx *>(b + covered);
          for (long start = 0; start < rest; ++start) {
            cplx v = c1 * a_c[start];
            if constexpr (Diag) v -= c2 * b_c[start];
            acc_c[start] += v;
          }
        }
      }
    }
  } // namespace

  void full_iw3ph_accumulate(const mc_weight_t sign, M3_M_full_t const &M, M3_GMG_t const &GMG, M3_G_t const &GM, M3_G_t const &MG, chi3_iw_v_t &M3,
                             const int bl1, const int bl2, const long bl2_size) noexcept {
    // Dispatch once per block pair, outside the mesh sweep: inside, N and diagonal are compile-time
    // so every accumulation loop has a known trip count and unrolls. A block size beyond max_block
    // maps to N == 0, the runtime-length path.
    const int n = (bl2_size <= max_block) ? static_cast<int>(bl2_size) : 0;
    poet::dispatch(
       [&]<int N, int Diag>() {
         constexpr bool diagonal = static_cast<bool>(Diag);

         // poet passes this lambda to the dispatch table by reference, so anything reached through a
         // capture stays a memory load the stores to acc can appear to alias. Copy it out once here.
         const mc_weight_t sgn = sign;

         auto const &M1 = M[bl1];
         // M2a is G*M*G: the one operand of this update that does not depend on the mesh point, so
         // its held register tiles below are loaded once for the whole sweep instead of per point.
         auto const &M2a     = GMG(bl2);
         auto const &GM1     = GM[bl1];
         auto const &MG2     = MG[bl2];
         auto const bl1_size = M1.target_shape()[0];
         auto &acc_bl        = M3(bl1, bl2);

         // acc(iw1, iw2) += sgn * M2a * M(iw1, iw2) [- sgn * GM1(iw1) * MG2(iw2) on the diagonal].
         // At N == 1 every block axis has length 1, so the loops below have nothing to vectorize and
         // would walk the mesh point by point. Both terms are separable in the two mesh indices, so
         // accumulate along the contiguous iw2 axis instead: one scaled add per acc row. Term 1's
         // coefficient is constant over the whole plane and M shares acc's mesh, so off the diagonal
         // it is a single scaled add over the flat plane.
         if constexpr (N == 1) {
           const long n1 = acc_bl.data().shape()[0], n2 = acc_bl.data().shape()[1];
           const cplx c1                        = sgn * M2a(0, 0);
           const auto *const TRIQS_RESTRICT m_d = reinterpret_cast<const double *>(M1.data().data());
           auto *const TRIQS_RESTRICT acc_d     = reinterpret_cast<double *>(acc_bl.data().data());
           if constexpr (diagonal) {
             const cplx *const TRIQS_RESTRICT gm1   = GM1.data().data();
             const auto *const TRIQS_RESTRICT mg2_d = reinterpret_cast<const double *>(MG2.data().data());
             const long row_len                     = 2 * n2;
             for (long i1 = 0; i1 < n1; ++i1) {
               const long row_start = i1 * row_len;
               fma_run<true, 0>(acc_d + row_start, m_d + row_start, mg2_d, c1, sgn * gm1[i1], row_len);
             }
           } else {
             fma_run<false, 0>(acc_d, m_d, m_d, c1, {}, 2 * n1 * n2);
           }
           return;
         }

         static_assert(std::decay_t<decltype(M2a.indexmap())>::is_stride_order_Fortran(),
                       "M2a must store its target transposed: it is read as M2a.data()[k*N + l]");
         const cplx *const TRIQS_RESTRICT M2a_flat = M2a.data();
         const auto *const TRIQS_RESTRICT M2a_d    = reinterpret_cast<const double *>(M2a_flat);

         // Run geometry. A length of 0 means the block size is only known at run time, which is what
         // fma_run's Len == 0 defers to its runtime argument.
         constexpr long ct_row_len   = 2 * long{N};    // doubles in one acc row
         constexpr long ct_plane_len = ct_row_len * N; // doubles in one acc(i,j) plane

         // Off the diagonal every (i,j) reads all of M2a, so hold it in vector registers across the loops
         // below: the plane as `tiles` widest vectors plus, when N is odd, one min_simd tail. The tile
         // indices must be compile-time (poet::static_for) or the array spills. At N == 8 this wants 16 of
         // the 32 zmm live and does spill, trading spill stores for not re-loading the operand on every
         // (i,j).
         constexpr bool hold       = N > 0 && !diagonal && worth_holding<N>;
         constexpr long tiles      = hold ? ct_plane_len / max_simd : 0;
         constexpr long tail_start = tiles * max_simd;
         constexpr long tail_len   = hold ? ct_plane_len % max_simd : 0;
         static_assert(tail_len == 0 || tail_len == min_simd, "the held plane needs a narrowest-batch tail to close an even run");
         std::array<vec<max_simd>, tiles> M2a_full;
         vec<min_simd> M2a_tail{};
         poet::static_for<tiles>([&]<auto Tile>() {
           constexpr auto start = Tile * max_simd;
           M2a_full[Tile]       = vec<max_simd>::load_unaligned(M2a_d + start);
         });
         if constexpr (tail_len) M2a_tail = vec<min_simd>::load_unaligned(M2a_d + tail_start);

         for (auto mp : acc_bl.mesh()) {
           auto [mp1, mp2] = mp;
           const auto iw1  = mp1.value();
           const auto iw2  = mp2.value();
           const auto M1a  = M1[closest_mesh_pt(iw1, iw2)];
           const auto M1b  = GM1[iw1];
           const auto M2b  = MG2[iw2];
           auto acc        = acc_bl[mp];

           static_assert(std::decay_t<decltype(M1b.indexmap())>::is_stride_order_Fortran(),
                         "M1b must store its target transposed: the column M1b(:,i) is read contiguously");
           static_assert(std::decay_t<decltype(acc.indexmap())>::is_stride_order_C(), "acc must be C-ordered: its (i,j) planes are walked linearly");
           assert(acc.indexmap().is_contiguous());

           // The acc(i,j,:,:) planes are contiguous and consecutive, so walk them with a
           // pointer. Re-deriving &acc(i,j,0,0) costs an nda 4-index offset plus view
           // refcounting, which at small N dominates the arithmetic itself.
           const long nrow            = (N > 0) ? long{N} : bl2_size; // rows in an acc(i,j) plane
           const long plane_len       = nrow * nrow;                  // complex in one plane
           cplx *TRIQS_RESTRICT plane = acc.data();
           using wide                 = vec<max_simd>;
           for (long i = 0; i < bl1_size; ++i) {
             const cplx *const TRIQS_RESTRICT M1b_col = &M1b(0, i);
             // The restrict-qualified bases stay outside the k loops below. Declared inside one,
             // each iteration opens a fresh restrict scope, so clang can no longer tell that the
             // store to acc leaves M1b(:,i) alone and re-loads plus re-swizzles it on every k.
             const auto *const TRIQS_RESTRICT M1b_d = reinterpret_cast<const double *>(M1b_col);
             for (long j = 0; j < bl1_size; ++j, plane += plane_len) {
               const cplx c1                    = M1a(j, i) * sgn;
               auto *const TRIQS_RESTRICT acc_d = reinterpret_cast<double *>(plane);

               if constexpr (hold) {
                 // plane(k,l) += c1 * M2a(l,k), the operand already sitting in registers.
                 const wide cr(c1.real()), ci(c1.imag());
                 poet::static_for<tiles>([&]<auto Tile>() {
                   constexpr auto start = Tile * max_simd;
                   cfma(cr, ci, M2a_full[Tile], wide::load_unaligned(acc_d + start)).store_unaligned(acc_d + start);
                 });
                 if constexpr (tail_len) {
                   using narrow = vec<min_simd>;
                   cfma(narrow(c1.real()), narrow(c1.imag()), M2a_tail, narrow::load_unaligned(acc_d + tail_start))
                      .store_unaligned(acc_d + tail_start);
                 }

               } else if constexpr (diagonal) {
                 // Row by row, so that M1b(:,i) is read once per row:
                 //   plane(k,:) += c1 * M2a(:,k) - c2(k) * M1b(:,i)
                 for (long k = 0; k < nrow; ++k) {
                   const long row_start = k * 2 * nrow;
                   fma_run<true, ct_row_len>(acc_d + row_start, M2a_d + row_start, M1b_d, c1, M2b(j, k) * sgn, 2 * nrow);
                 }

               } else {
                 // The plane is contiguous, so its rows are one run: plane += c1 * M2a.
                 fma_run<false, ct_plane_len>(acc_d, M2a_d, M1b_d, c1, {}, 2 * plane_len);
               }
             }
           }
         }
       },
       poet::dispatch_param<poet::inclusive_range<0, max_block>>{n}, poet::dispatch_param<poet::inclusive_range<0, 1>>{bl1 == bl2});
  }

  M3ph_iw_full::M3ph_iw_full(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_),
       qmc_config(qmc_config_),
       buf_arrarr(params_.n_blocks()),
       buf_arrarr_GM(params_.n_blocks()),
       buf_arrarr_MG(params_.n_blocks()),
       G0_tau(std::move(G0_tau_)) {

    // Construct full fermionic Matsubara mesh
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M3};
    mesh::prod iw2_mesh{iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M3ph_iw_nfft_full = make_block2_gf(iw2_mesh, params.gf_struct);
    M3ph_iw_.rebind(results->M3ph_iw_nfft_full.value());
    M3ph_iw_() = 0;

    // Initialize M on full 2D uniform grid (type1 Rank=2 NFFT)
    M = block_gf{iw2_mesh, params.gf_struct};

    // Initialize GM, MG on uniform imfreq mesh (for type1 NFFT)
    GM = block_gf{iw_mesh, params.gf_struct};
    MG = block_gf{iw_mesh, params.gf_struct};

    auto init_target_func = [&](int bl) {
      int bl_size = GM[bl].target_shape()[0];
      return array<dcomplex, 2>(bl_size, bl_size);
    };
    GMG = array_adapter{make_shape(params.n_blocks()), init_target_func};

    // Create nfft buffers: type1 Rank=2 for M (uniform 2D grid), type1 Rank=1 for GM/MG
    for (auto bl : range(params.n_blocks())) {
      buf_arrarr(bl) =
         array_adapter{M[bl].target_shape(), [&](int i, int j) {
                         return nfft::buffer_t{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
                       }};
      buf_arrarr_GM(bl) =
         array_adapter{GM[bl].target_shape(), [&](int i, int j) {
                         return nfft::buffer_t{slice_target_to_scalar(GM[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
                       }};
      buf_arrarr_MG(bl) =
         array_adapter{MG[bl].target_shape(), [&](int i, int j) {
                         return nfft::buffer_t{slice_target_to_scalar(MG[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
                       }};
    }
  }

  void M3ph_iw_full::accumulate(mc_weight_t sign) {
    Z += sign;

    // Reset all accumulators
    M()  = 0;
    GM() = 0;
    MG() = 0;

    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = GM[bl].target_shape()[0];
      long k      = det.size();
      if (k == 0) {
        GMG(bl)() = 0;
        continue;
      }

      auto Ginv = det.inverse_matrix();

      // Single-orbital fast path: fused (i, j) loop over M, GMG, GM, MG avoids
      // GEMV overhead that dominates at bl_size=1 where GEMMs do not amortize.
      if (bl_size == 1) {
        auto arr_GM = nda::zeros<dcomplex>(k);
        auto arr_MG = nda::zeros<dcomplex>(k);
        GMG(bl)()   = 0;
        for (long i = 0; i < k; ++i) {
          auto tau_i = double(det.get_x(i).tau);
          auto G0_i  = G0_tau[bl][closest_mesh_pt(tau_i)](0, 0);
          for (long j = 0; j < k; ++j) {
            auto tau_j = double(det.get_y(j).tau);
            auto G0_j  = -G0_tau[bl][closest_mesh_pt(params.beta - tau_j)](0, 0);
            auto g_ji  = Ginv(j, i);
            // Coordinates in (beta - tau_i, tau_j) order so the NFFT lands M already transposed:
            // the kernel wants M(iw1, iw2) and would otherwise have to read it strided.
            buf_arrarr(bl)(0, 0).push_back({params.beta - tau_i, tau_j}, -g_ji);
            GMG(bl)(0, 0) += G0_j * g_ji * G0_i;
            arr_GM(i) += -G0_j * g_ji;
            arr_MG(j) += g_ji * G0_i;
          }
        }
        for (long p = 0; p < k; ++p) {
          buf_arrarr_GM(bl)(0, 0).push_back({params.beta - double(det.get_x(p).tau)}, arr_GM(p));
          buf_arrarr_MG(bl)(0, 0).push_back({double(det.get_y(p).tau)}, arr_MG(p));
        }
        for (auto &buf : buf_arrarr(bl)) buf.flush();
        for (auto &buf : buf_arrarr_GM(bl)) buf.flush();
        for (auto &buf : buf_arrarr_MG(bl)) buf.flush();
        continue;
      }

      // M on full 2D grid via type1 NFFT
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (long j = 0; j < k; ++j) {
          auto &[tau_j, u_j, _, _, _] = det.get_y(j);
          // See the bl_size == 1 branch: transposed coordinate order, so M arrives as M(iw1, iw2).
          buf_arrarr(bl)(u_j, u_i).push_back({params.beta - double(tau_i), double(tau_j)}, -Ginv(j, i));
        }
      }
      for (auto &buf : buf_arrarr(bl)) buf.flush();

      // Build X(a, i) = G0(tau_i)(u_i, a) and Y(b, j) = -G0(beta - tau_j)(b, u_j)
      auto X = nda::matrix<dcomplex>(bl_size, k);
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        auto G0_i                   = G0_tau[bl][closest_mesh_pt(double(tau_i))];
        for (int a = 0; a < bl_size; ++a) X(a, i) = G0_i(u_i, a);
      }

      auto Y = nda::matrix<dcomplex>(bl_size, k);
      for (long j = 0; j < k; ++j) {
        auto &[tau_j, u_j, _, _, _] = det.get_y(j);
        auto G0_j                   = G0_tau[bl][closest_mesh_pt(params.beta - tau_j)];
        for (int b = 0; b < bl_size; ++b) Y(b, j) = -G0_j(b, u_j);
      }

      // W(b, i) = sum_j Y(b, j) Ginv(j, i): the tau-space GM = G_left * M (G_left carries the minus in Y)
      auto W = Y * Ginv;

      // GMG(a, b) = sum_i W(a, i) X(b, i) = (G_left * M * G_right)(a, b), with a the G_left and b the
      // G_right orbital. The accumulation kernel reads GMG(l, k), so this (W * X^T) is the required order.
      GMG(bl) = W * nda::transpose(X);

      // GM: scatter to NFFT buffers. The kernel consumes the same GM as M3pp (built from +G0(beta - tau_j)),
      // which is -W here since Y carries the extra minus. One value per orbital pair, no bl_size factor.
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (int m : range(bl_size)) buf_arrarr_GM(bl)(m, u_i).push_back({params.beta - double(tau_i)}, -W(m, i));
      }
      for (auto &buf : buf_arrarr_GM(bl)) buf.flush();

      // MG: V(j, n) = sum_i Ginv(j, i) X(n, i) = (M * G_right)(j, n). Scatter to row u_j, keeping the free
      // orbital n (the second G0 index) intact -- it must not be summed over.
      auto V = Ginv * nda::transpose(X);
      for (long j = 0; j < k; ++j) {
        auto u_j = det.get_y(j).u;
        for (int n : range(bl_size)) buf_arrarr_MG(bl)(u_j, n).push_back({double(det.get_y(j).tau)}, V(j, n));
      }
      for (auto &buf : buf_arrarr_MG(bl)) buf.flush();
    }

    // Accumulation kernel
    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) full_iw3ph_accumulate(sign, M, GMG, GM, MG, M3ph_iw_, bl1, bl2, GMG(bl2).shape()[0]);
  }

  void M3ph_iw_full::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3ph_iw_ = mpi::all_reduce(M3ph_iw_, comm);
    M3ph_iw_ = M3ph_iw_ / Z;
  }

} // namespace triqs_ctint::measures
