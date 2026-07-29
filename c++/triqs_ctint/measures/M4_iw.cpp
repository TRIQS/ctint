// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M4_iw.hpp"
#include "./iw_simd.hpp"

#include <poet/poet.hpp>

#include <array>
#include <cassert>

namespace triqs_ctint::measures {

  void iw4_accumulate(const mc_weight_t sign, M4_M_t const &M, chi4_iw_v_t &M4, const int bl1, const int bl2, const long bl2_size) noexcept {
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

         // All three components of M4's mesh are the same iw mesh; bind it once and walk it thrice.
         auto &[iw_mesh, _, _] = M4(0, 0).mesh();
         auto const &M1        = M[bl1];
         auto const &M2        = M[bl2];
         auto const bl1_size   = M1.target_shape()[0];
         auto &acc_bl          = M4(bl1, bl2);

         for (auto const &iw1 : iw_mesh) {
           for (auto const &iw2 : iw_mesh) {
             for (auto const &iw3 : iw_mesh) {
               const auto iw4 = iw1 + iw3 - iw2;
               const auto M1a = M1[iw2.value(), iw1];
               const auto M2a = M2[iw4, iw3];
               const auto M1b = M1[iw4, iw1];
               const auto M2b = M2[iw2.value(), iw3];
               auto acc       = acc_bl[iw1, iw2, iw3];

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

               static_assert(std::decay_t<decltype(M1b.indexmap())>::is_stride_order_Fortran(),
                             "M1b must store its target transposed: the column M1b(:,i) is read contiguously");
               static_assert(std::decay_t<decltype(acc.indexmap())>::is_stride_order_C(),
                             "acc must be C-ordered: its (i,j) planes are walked linearly");
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
                         // 2*N contiguous doubles, covered widest-batch-first.
                         const long base = long{k} * 2 * long{N};
                         simd_cover<2 * long{N}>([&](auto w, const long off) {
                           using batch = typename decltype(w)::type;
                           auto av     = batch::load_unaligned(acc_d + base + off);
                           av += cmul(batch(c1.real()), batch(c1.imag()), batch::load_unaligned(M2a_d + base + off))
                              - cmul(batch(c2.real()), batch(c2.imag()), batch::load_unaligned(M1b_d + off));
                           av.store_unaligned(acc_d + base + off);
                         });
                       } else
                         for (int l = 0; l < N; ++l) plane[k * N + l] += c1 * M2a_flat[k * N + l] - c2 * M1b_col[l];
                     }
                   }
                 }
               }
             }
           }
         }
       },
       poet::dispatch_param<poet::inclusive_range<0, max_block>>{n}, poet::dispatch_param<poet::inclusive_range<0, 1>>{bl1 == bl2});
  }

  M4_iw::M4_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M4};
    mesh::prod<imfreq, imfreq, imfreq> M4_iw_mesh{iw_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M4_iw = make_block2_gf(M4_iw_mesh, params.gf_struct);
    M4_iw_.rebind(results->M4_iw.value());
    M4_iw_() = 0;

    // Construct Matsubara mesh for temporary Matrix
    mesh::imfreq iw_mesh_large{params.beta, Fermion, 3 * params.n_iw_M4};
    mesh::prod<imfreq, imfreq> M_mesh{iw_mesh_large, iw_mesh};

    // Initialize intermediate scattering matrix
    M = block_gf{M_mesh, params.gf_struct};

    // Create nfft buffers
    for (int bl : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft::buffer_t<2>{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
      };
      buf_arrarr(bl) = array_adapter{M[bl].target_shape(), init_target_func};
    }
  }

  void M4_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    M() = 0;
    for (int bl : range(params.n_blocks()))
      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[bl],
               [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) { // Care for negative frequency in c transform (for M-objects)
                 buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({double(cdag_j.tau), params.beta - double(c_i.tau)}, -Ginv_ji);
               });
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) iw4_accumulate(sign, M, M4_iw_, bl1, bl2, M[bl2].target_shape()[0]);
  }

  void M4_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z      = mpi::all_reduce(Z, comm);
    M4_iw_ = mpi::all_reduce(M4_iw_, comm);
    M4_iw_ = M4_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
