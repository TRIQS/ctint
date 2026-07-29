// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3ph_iw.hpp"
#include "./iw_simd.hpp"

#include <poet/poet.hpp>

#include <array>
#include <cassert>

namespace triqs_ctint::measures {

  void dlr2d_iw3ph_accumulate(const mc_weight_t sign, M3_M_t const &M, M3_GMG_t const &GMG, M3_G_t const &GM, M3_G_t const &MG, chi3_dlr2d_iw_v_t &M3,
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

         static_assert(std::decay_t<decltype(M2a.indexmap())>::is_stride_order_Fortran(),
                       "M2a must store its target transposed: it is read as M2a.data()[k*N + l]");
         const cplx *const TRIQS_RESTRICT M2a_flat = M2a.data();
         const auto *const TRIQS_RESTRICT M2a_d    = reinterpret_cast<const double *>(M2a_flat);

         // Off the diagonal every (i,j) reads all N*N of M2a, so hold it in vector registers
         // across the loops below: 2*N*N doubles as `held_full` widest tiles plus, when N is
         // odd, one min_simd tail. The tile indices must be compile-time (poet::static_for) or
         // the array spills. At N == 8 this wants 16 of the 32 zmm live and does spill,
         // trading spill stores for not re-loading the operand on every (i,j).
         constexpr bool hold      = N > 0 && !diagonal && prefer_simd<N>;
         constexpr long held_full = hold ? 2 * long{N} * N / max_simd : 0;
         constexpr long held_tail = hold ? 2 * long{N} * N % max_simd : 0;
         static_assert(held_tail == 0 || held_tail == min_simd);
         std::array<vec<max_simd>, held_full> M2a_full;
         vec<min_simd> M2a_tail{};
         if constexpr (hold) {
           poet::static_for<held_full>([&]<auto Tile>() { M2a_full[Tile] = vec<max_simd>::load_unaligned(M2a_d + Tile * max_simd); });
           if constexpr (held_tail) M2a_tail = vec<min_simd>::load_unaligned(M2a_d + held_full * max_simd);
         }

         for (auto mp : acc_bl.mesh()) {
           auto [iw1, iw2] = mp.value();
           const auto M1a  = M1[mp];
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
           const long plane_len       = (N > 0) ? long{N} * N : bl2_size * bl2_size;
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

               if constexpr (N == 0) {
                 // Runtime block size: whole widest vectors, then a scalar tail. max_simd is a
                 // power of two, so the mask truncates to the last whole batch.
                 if constexpr (!diagonal) {
                   const long lanes = (2 * plane_len) & -max_simd;
                   const wide cr(c1.real()), ci(c1.imag());
                   for (long off = 0; off < lanes; off += max_simd)
                     cfma(cr, ci, wide::load_unaligned(M2a_d + off), wide::load_unaligned(acc_d + off)).store_unaligned(acc_d + off);
                   for (long l = lanes / 2; l < plane_len; ++l) plane[l] += c1 * M2a_flat[l];
                 } else {
                   const long lanes = (2 * bl2_size) & -max_simd;
                   const wide c1r(c1.real()), c1i(c1.imag());
                   for (long k = 0; k < bl2_size; ++k) {
                     const cplx c2 = M2b(j, k) * sgn;
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
                 // plane[0 : N*N) += c1 * M2a, the operand already sitting in registers.
                 if constexpr (prefer_simd<N>) {
                   const wide cr(c1.real()), ci(c1.imag());
                   poet::static_for<held_full>([&]<auto Tile>() {
                     cfma(cr, ci, M2a_full[Tile], wide::load_unaligned(acc_d + Tile * max_simd)).store_unaligned(acc_d + Tile * max_simd);
                   });
                   if constexpr (held_tail) {
                     using narrow       = vec<min_simd>;
                     constexpr long off = held_full * max_simd;
                     cfma(narrow(c1.real()), narrow(c1.imag()), M2a_tail, narrow::load_unaligned(acc_d + off)).store_unaligned(acc_d + off);
                   }
                 } else
                   for (long l = 0; l < long{N} * N; ++l) plane[l] += c1 * M2a_flat[l];

               } else {
                 for (int k = 0; k < N; ++k) {
                   const cplx c2 = M2b(j, k) * sgn;
                   // plane[k*N : k*N + N) += c1 * M2a(:,k) - c2 * M1b(:,i), in one pass over acc.
                   if constexpr (prefer_simd<N>) {
                     // 2*N contiguous doubles: as many widest ops as fit, then the even
                     // remainder as at most one 4-wide plus one 2-wide op, each writing exactly
                     // the bytes it covers. A single masked widest store would cover the
                     // remainder in far fewer instructions but writes a subset of a widest range
                     // that the next k's load overlaps, blocking store-to-load forwarding.
                     static_assert(min_simd == 2 && max_simd <= 8, "the remainder cover below enumerates the cases for min_simd == 2, max_simd <= 8");
                     constexpr long len = 2 * long{N}, rem = len % max_simd;
                     const long base = long{k} * len;
                     auto op         = [&](auto w, const long off) {
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
       },
       poet::dispatch_param<poet::inclusive_range<0, max_block>>{n}, poet::dispatch_param<poet::inclusive_range<0, 1>>{bl1 == bl2});
  }

  M3ph_iw::M3ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), buf_arrarr_GM(params_.n_blocks()), buf_arrarr_MG(params_.n_blocks()), G0_tau(std::move(G0_tau_)) {

    // Construct DLR2D Matsubara mesh
    mesh::dlr2d_imfreq M3ph_iw_mesh{params.beta, params.dlr_wmax, params.dlr_eps, mesh::PH, params.dlr2d_compress_grid};

    // Init measurement container and capture view
    results->M3ph_iw_nfft = make_block2_gf(M3ph_iw_mesh, params.gf_struct);
    M3ph_iw_.rebind(results->M3ph_iw_nfft.value());
    M3ph_iw_() = 0;

    // Initialize M on DLR2D mesh
    M = block_gf{M3ph_iw_mesh, params.gf_struct};

    // Build target frequencies for matrix_buffer_t (PH convention: {iw2, iw1})
    std::vector<std::array<mesh::matsubara_freq, 2>> target_mf_2d;
    target_mf_2d.reserve(M3ph_iw_mesh.size());
    for (auto mp : M3ph_iw_mesh) {
      auto [iw1, iw2] = mp.value();
      target_mf_2d.push_back({iw2, iw1});
    }

    // Create matrix buffers for M (factored product-grid DFT)
    M_bufs.reserve(params.n_blocks());
    for (auto bl : range(params.n_blocks())) {
      int bl_size = params.gf_struct[bl].second;
      M_bufs.emplace_back(M[bl].data(), target_mf_2d, bl_size);
    }

    // Initialize GM, MG on uniform imfreq mesh (for type1 NFFT)
    mesh::imfreq iw_mesh{params.beta, Fermion, M3ph_iw_mesh.max_n() + 1};
    GM = block_gf{iw_mesh, params.gf_struct};
    MG = block_gf{iw_mesh, params.gf_struct};

    auto init_target_func = [&](int bl) {
      int bl_size = GM[bl].target_shape()[0];
      return array<dcomplex, 2>(bl_size, bl_size);
    };
    GMG = array_adapter{make_shape(params.n_blocks()), init_target_func};

    // Create nfft buffers: type1 for GM/MG (uniform grid)
    for (auto bl : range(params.n_blocks())) {
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

  void M3ph_iw::accumulate(mc_weight_t sign) {
    Z += sign;

    // Reset all accumulators
    M()  = 0;
    GM() = 0;
    MG() = 0;

    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = M[bl].target_shape()[0];
      long k      = det.size();
      if (k == 0) {
        GMG(bl)() = 0;
        continue;
      }

      auto Ginv = det.inverse_matrix();

      // M on DLR2D mesh via factored product-grid DFT
      std::vector<double> x(k), y(k);
      std::vector<int> u_x(k), u_y(k);
      for (long j = 0; j < k; ++j) {
        x[j]   = double(det.get_y(j).tau);
        u_x[j] = det.get_y(j).u;
      }
      for (long i = 0; i < k; ++i) {
        y[i]   = params.beta - double(det.get_x(i).tau);
        u_y[i] = det.get_x(i).u;
      }
      M_bufs[bl].push_product(x.data(), u_x.data(), k, y.data(), u_y.data(), k, nda::matrix<dcomplex>(-Ginv));

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
      for (int bl2 : range(params.n_blocks())) dlr2d_iw3ph_accumulate(sign, M, GMG, GM, MG, M3ph_iw_, bl1, bl2, GMG(bl2).shape()[0]);
  }

  void M3ph_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3ph_iw_ = mpi::all_reduce(M3ph_iw_, comm);
    M3ph_iw_ = M3ph_iw_ / Z;
  }

} // namespace triqs_ctint::measures
