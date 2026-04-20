// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3ph_iw.hpp"
#include "./iw_accumulate.hpp"

namespace triqs_ctint::measures {

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
    for (auto &gmg : GMG) gmg() = 0;
    M()  = 0;
    GM() = 0;
    MG() = 0;

    double beta = params.beta;

    // Phase 1: M construction via factored product-grid DFT
    for (int bl : range(params.n_blocks())) {
      auto &det = qmc_config.dets[bl];
      long k    = det.size();
      if (k == 0) continue;

      auto Ginv = det.inverse_matrix();

      // Compute M on DLR2D mesh via factored product-grid DFT
      std::vector<double> x(k), y(k);
      std::vector<int> u_x(k), u_y(k);
      for (long j = 0; j < k; ++j) {
        x[j]   = double(det.get_y(j).tau);
        u_x[j] = det.get_y(j).u;
      }
      for (long i = 0; i < k; ++i) {
        y[i]   = beta - double(det.get_x(i).tau);
        u_y[i] = det.get_x(i).u;
      }
      auto neg_Ginv = nda::matrix<dcomplex>(-Ginv);
      M_bufs[bl].push_product(x.data(), u_x.data(), k, y.data(), u_y.data(), k, neg_Ginv);
    }

    // Phase 2: GM/MG/GMG construction (BLAS GEMM + NFFT scatter + flush)
    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = GM[bl].target_shape()[0];
      long k      = det.size();
      if (k == 0) continue;

      auto Ginv = det.inverse_matrix();

      auto X = nda::matrix<dcomplex>(bl_size, k);
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        auto G0_i                   = G0_tau[bl][closest_mesh_pt(double(tau_i))];
        for (int a = 0; a < bl_size; ++a) X(a, i) = G0_i(u_i, a);
      }

      auto Y = nda::matrix<dcomplex>(bl_size, k);
      for (long j = 0; j < k; ++j) {
        auto &[tau_j, u_j, _, _, _] = det.get_y(j);
        auto G0_j                   = G0_tau[bl][closest_mesh_pt(beta - tau_j)];
        for (int b = 0; b < bl_size; ++b) Y(b, j) = -G0_j(b, u_j);
      }

      auto W  = Y * Ginv;
      GMG(bl) = X * nda::transpose(W);

      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (int m : range(bl_size))
          buf_arrarr_GM(bl)(m, u_i).push_back({beta - double(tau_i)}, dcomplex(-bl_size) * W(m, i));
      }

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
          for (int m : range(bl_size))
            buf_arrarr_MG(bl)(m, n).push_back({tau_j}, MG_temp(n, j));
      }
    }

    for (auto &buf_arr : buf_arrarr_GM)
      for (auto &buf : buf_arr) buf.flush();
    for (auto &buf_arr : buf_arrarr_MG)
      for (auto &buf : buf_arr) buf.flush();

    // Phase 3: Accumulation kernel
    for (auto bl1 : range(params.n_blocks()))
      for (auto bl2 : range(params.n_blocks())) {
        auto const bl2_size = GMG(bl2).shape()[0];
        simd::dlr2d_iw3ph_accumulate(sign, M, GMG, GM, MG, M3ph_iw_, bl1, bl2, bl2_size);
      }
  }

  void M3ph_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3ph_iw_ = mpi::all_reduce(M3ph_iw_, comm);
    M3ph_iw_ = M3ph_iw_ / Z;
  }

} // namespace triqs_ctint::measures
