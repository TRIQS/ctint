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
    for (auto &gmg : GMG) gmg() = 0;
    M()  = 0;
    GM() = 0;
    MG() = 0;

    double beta = params.beta;

    // Compute M via type1 NFFT, and GMG, GM, MG via BLAS GEMM
    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = GM[bl].target_shape()[0];
      long k      = det.size();
      if (k == 0) continue;

      // Get full inverse matrix
      auto Ginv = det.inverse_matrix(); // k x k

      // Build X(a, i) = G0(tau_i)(u_i, a) -- one G0 lookup per i
      auto X = nda::matrix<dcomplex>(bl_size, k);
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        auto G0_i                   = G0_tau[bl][closest_mesh_pt(double(tau_i))];
        for (int a = 0; a < bl_size; ++a) X(a, i) = G0_i(u_i, a);
      }

      // Build Y(b, j) = -G0(beta - tau_j)(b, u_j) -- one G0 lookup per j
      auto Y = nda::matrix<dcomplex>(bl_size, k);
      for (long j = 0; j < k; ++j) {
        auto &[tau_j, u_j, _, _, _] = det.get_y(j);
        auto G0_j                   = G0_tau[bl][closest_mesh_pt(beta - tau_j)];
        for (int b = 0; b < bl_size; ++b) Y(b, j) = -G0_j(b, u_j);
      }

      // Compute M on full 2D grid via type1 NFFT
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (long j = 0; j < k; ++j) {
          auto &[tau_j, u_j, _, _, _] = det.get_y(j);
          buf_arrarr(bl)(u_j, u_i).push_back({double(tau_j), beta - double(tau_i)}, -Ginv(j, i));
        }
      }

      // W = Y * Ginv via BLAS GEMM: (bl_size, k) * (k, k) -> (bl_size, k)
      auto W = Y * Ginv;

      // GMG = X * W^T via BLAS GEMM: (bl_size, k) * (k, bl_size) -> (bl_size, bl_size)
      GMG(bl) = X * nda::transpose(W);

      // GM: scatter W into NFFT buffers (factor -bl_size from redundant abar_u loop in original)
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (int m : range(bl_size)) buf_arrarr_GM(bl)(m, u_i).push_back({beta - double(tau_i)}, dcomplex(-bl_size) * W(m, i));
      }

      // MG: scatter XsumScatter * Ginv^T into NFFT buffers
      auto XsumScatter = nda::matrix<dcomplex>(bl_size, k);
      XsumScatter()    = 0;
      for (long i = 0; i < k; ++i) {
        auto u_i    = det.get_x(i).u;
        dcomplex xs = 0;
        for (int a = 0; a < bl_size; ++a) xs += X(a, i);
        XsumScatter(u_i, i) = xs;
      }
      auto MG_temp = XsumScatter * nda::transpose(Ginv);

      for (long j = 0; j < k; ++j) {
        auto tau_j = double(det.get_y(j).tau);
        for (int n : range(bl_size))
          for (int m : range(bl_size)) buf_arrarr_MG(bl)(m, n).push_back({tau_j}, MG_temp(n, j));
      }
    }

    // Flush remaining points from all NFFT buffers
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush();
    for (auto &buf_arr : buf_arrarr_GM)
      for (auto &buf : buf_arr) buf.flush();
    for (auto &buf_arr : buf_arrarr_MG)
      for (auto &buf : buf_arr) buf.flush();

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
