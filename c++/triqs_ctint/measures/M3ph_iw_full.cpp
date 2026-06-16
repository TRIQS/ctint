// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3ph_iw_full.hpp"
#include "./iw_accumulate.hpp"

namespace triqs_ctint::measures {

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
            buf_arrarr(bl)(0, 0).push_back({tau_j, params.beta - tau_i}, -g_ji);
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
          buf_arrarr(bl)(u_j, u_i).push_back({double(tau_j), params.beta - double(tau_i)}, -Ginv(j, i));
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
    for (auto bl1 : range(params.n_blocks()))
      for (auto bl2 : range(params.n_blocks())) {
        auto const bl2_size = GMG(bl2).shape()[0];
        simd::full_iw3ph_accumulate(sign, M, GMG, GM, MG, M3ph_iw_, bl1, bl2, bl2_size);
      }
  }

  void M3ph_iw_full::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3ph_iw_ = mpi::all_reduce(M3ph_iw_, comm);
    M3ph_iw_ = M3ph_iw_ / Z;
  }

} // namespace triqs_ctint::measures
