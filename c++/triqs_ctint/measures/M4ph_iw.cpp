// Copyright (c) 2023--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M4ph_iw.hpp"
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

  void iw4ph_accumulate(const mc_weight_t sign, M4_M_t const &M, chi4_iw_v_t &M4ph, const int bl1, const int bl2, const long bl2_size) noexcept {
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

         auto const &[iW_mesh, iw_mesh, _] = M4ph(0, 0).mesh();
         auto const &M1                    = M[bl1];
         auto const &M2                    = M[bl2];
         auto const bl1_size               = M1.target_shape()[0];
         auto &acc_bl                      = M4ph(bl1, bl2);

         for (auto iW : iW_mesh) {
           for (auto iw : iw_mesh) {
             // iW + iw is fixed for the whole iwp sweep, and so is the M1a it selects.
             const auto iW_iw = iW + iw;
             const auto M1a   = M1[iW_iw, iw.value()];
             for (auto iwp : iw_mesh) {
               const auto iW_iwp = iW + iwp;
               const auto M2a    = M2[iwp.value(), iW_iwp];
               const auto M1b    = M1[iwp.value(), iw.value()];
               const auto M2b    = M2[iW_iw, iW_iwp];
               auto acc          = acc_bl[iW, iw, iwp];

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

               static_assert(std::decay_t<decltype(M1b.indexmap())>::is_stride_order_Fortran(),
                             "M1b must store its target transposed: the column M1b(:,i) is read contiguously");
               static_assert(std::decay_t<decltype(acc.indexmap())>::is_stride_order_C(),
                             "acc must be C-ordered: its (i,j) planes are walked linearly");
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
           }
         }
       },
       poet::dispatch_param<poet::inclusive_range<0, max_block>>{n}, poet::dispatch_param<poet::inclusive_range<0, 1>>{bl1 == bl2});
  }

  M4ph_iw::M4ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    mesh::imfreq iW_mesh{params.beta, Boson, params.n_iW_M4};
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M4};
    mesh::prod<imfreq, imfreq, imfreq> M4ph_iw_mesh{iW_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M4ph_iw = make_block2_gf(M4ph_iw_mesh, params.gf_struct);
    M4ph_iw_.rebind(results->M4ph_iw.value());
    M4ph_iw_() = 0;

    // Construct Matsubara mesh for temporary Matrix
    mesh::imfreq iw_mesh_large{params.beta, Fermion, params.n_iW_M4 + params.n_iw_M4 + 1};
    mesh::prod<imfreq, imfreq> M_mesh{iw_mesh_large, iw_mesh_large};

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

  void M4ph_iw::accumulate(mc_weight_t sign) {
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
      for (int bl2 : range(params.n_blocks())) iw4ph_accumulate(sign, M, M4ph_iw_, bl1, bl2, M[bl2].target_shape()[0]);
  }

  void M4ph_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M4ph_iw_ = mpi::all_reduce(M4ph_iw_, comm);
    M4ph_iw_ = M4ph_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
