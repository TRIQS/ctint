// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3ph_iw.hpp"

namespace triqs_ctint::measures {

  M3ph_iw::M3ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_),
       qmc_config(qmc_config_),
       buf_arrarr(params_.n_blocks()),
       buf_arrarr_GM(params_.n_blocks()),
       buf_arrarr_MG(params_.n_blocks()),
       G0_tau(std::move(G0_tau_)) {

    // Construct DLR2D Matsubara mesh
    mesh::dlr2d_imfreq M3ph_iw_mesh{params.beta, params.dlr_wmax_M3, params.dlr_eps_M3, mesh::PH};

    // Init measurement container and capture view
    results->M3ph_iw_nfft = make_block2_gf(M3ph_iw_mesh, params.gf_struct);
    M3ph_iw_.rebind(results->M3ph_iw_nfft.value());
    M3ph_iw_() = 0;

    // Initialize intermediate scattering matrices on regular imfreq mesh (for type1 NFFT)
    mesh::imfreq iw_mesh{params.beta, Fermion, M3ph_iw_mesh.max_n() + 1};
    M  = block_gf{mesh::prod<imfreq, imfreq>{iw_mesh, iw_mesh}, params.gf_struct};
    GM = block_gf{iw_mesh, params.gf_struct};
    MG = block_gf{iw_mesh, params.gf_struct};

    auto init_target_func = [&](int bl) {
      int bl_size = GM[bl].target_shape()[0];
      return array<dcomplex, 2>(bl_size, bl_size);
    };
    GMG = array_adapter{make_shape(params.n_blocks()), init_target_func};

    // Create type1 nfft buffers that write to the block_gf data
    for (auto bl : range(params.n_blocks())) {
      buf_arrarr(bl) = array_adapter{M[bl].target_shape(), [&](int i, int j) {
        return nfft_buf_t{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
      }};
      buf_arrarr_GM(bl) = array_adapter{GM[bl].target_shape(), [&](int i, int j) {
        return nfft_buf_t{slice_target_to_scalar(GM[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
      }};
      buf_arrarr_MG(bl) = array_adapter{MG[bl].target_shape(), [&](int i, int j) {
        return nfft_buf_t{slice_target_to_scalar(MG[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
      }};
    }
  }

  void M3ph_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Reset intermediate scattering matrices
    for (auto &gmg : GMG) { gmg() = 0; }
    M()  = 0;
    GM() = 0;
    MG() = 0;

    double beta = params.beta;

    // Fill intermediate scattering matrices using explicit row-column loops
    // to accumulate values before pushing to NFFT buffers:
    // - GM: for fixed c_i (row), sum over all cdag_j (columns) before pushing once
    // - MG: for fixed cdag_j (col), sum over all c_i (rows) before pushing once
    for (int bl : range(params.n_blocks())) {
      int bl_size = params.gf_struct[bl].second;
      auto &det   = qmc_config.dets[bl];
      long k_bl   = det.size();

      nda::array<dcomplex, 1> gm_acc(bl_size);       // GM accumulator per b_u for current row
      nda::array<dcomplex, 2> mg_acc(k_bl, bl_size); // MG accumulator per (col, orbital)
      mg_acc = 0;

      for (long row = 0; row < k_bl; ++row) {
        auto const &c_i = det.get_x(row);
        auto tau_i      = double(c_i.tau);
        auto G0_tau_i   = G0_tau[bl][closest_mesh_pt(tau_i)];

        // Pre-compute sum_a G0(tau_i)(c_i.u, a) for MG
        dcomplex mg_g0_sum = 0;
        for (int a : range(bl_size)) mg_g0_sum += G0_tau_i(c_i.u, a);

        gm_acc = 0;

        for (long col = 0; col < k_bl; ++col) {
          auto const &cdag_j = det.get_y(col);
          auto Ginv_ji       = det.inverse_matrix(col, row);
          auto tau_j         = double(cdag_j.tau);

          // Fill M, Note: Minus sign from the shift of -tau_i
          buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({tau_j, beta - tau_i}, -Ginv_ji);

          // Accumulate GM and fill GMG
          for (int b_u : range(bl_size)) {
            // Note: Minus sign from the shift of -tau_j
            auto G0_bj = -G0_tau[bl][closest_mesh_pt(beta - tau_j)](b_u, cdag_j.u);
            gm_acc(b_u) += -G0_bj * Ginv_ji;
            for (int abar_u : range(bl_size)) {
              auto G0_ia = G0_tau_i(c_i.u, abar_u);
              GMG(bl)(abar_u, b_u) += G0_bj * Ginv_ji * G0_ia;
            }
          }

          // Accumulate MG: sum over rows for this column and orbital
          mg_acc(col, c_i.u) += Ginv_ji * mg_g0_sum;
        }

        // Push accumulated GM once per (b_u, c_i.u) — replaces k_bl * bl_size pushes
        // Note: Minus sign from the shift of -tau_i; bl_size factor from abar_u summation
        for (int b_u : range(bl_size)) buf_arrarr_GM(bl)(b_u, c_i.u).push_back({beta - tau_i}, double(bl_size) * gm_acc(b_u));
      }

      // Push accumulated MG once per (col, b_u, orbital) — replaces k_bl * bl_size pushes per column
      for (long col = 0; col < k_bl; ++col) {
        auto tau_j = double(det.get_y(col).tau);
        for (int u_i : range(bl_size))
          for (int b_u : range(bl_size)) buf_arrarr_MG(bl)(b_u, u_i).push_back({tau_j}, mg_acc(col, u_i));
      }
    }

    // Flush remaining points from all buffers
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush();
    for (auto &buf_arr : buf_arrarr_GM)
      for (auto &buf : buf_arr) buf.flush();
    for (auto &buf_arr : buf_arrarr_MG)
      for (auto &buf : buf_arr) buf.flush();

    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        int bl1_size     = M[bl1].target_shape()[0];
        int bl2_size     = M[bl2].target_shape()[0];
        auto const &M1   = M[bl1];
        auto const &GMG2 = GMG(bl2);
        auto const &GM1  = GM[bl1];
        auto const &MG2  = MG[bl2];
        auto &M3ph_iw    = M3ph_iw_(bl1, bl2);

        // Single loop over DLR2D mesh points
        for (auto mp : M3ph_iw.mesh()) {
          auto [iw1, iw2] = mp.value(); // matsubara_freq pair
          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size)) {
                  // Note: PH channel uses transposed access M1[iw2, iw1]
                  M3ph_iw[mp](i, j, k, l) += sign * M1[iw2, iw1](j, i) * GMG2(l, k);
                  if (bl1 == bl2) { M3ph_iw[mp](i, j, k, l) -= sign * GM1[iw1](l, i) * MG2[iw2](j, k); }
                }
        }
      }
  }

  void M3ph_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3ph_iw_ = mpi::all_reduce(M3ph_iw_, comm);
    M3ph_iw_ = M3ph_iw_ / Z;
  }

} // namespace triqs_ctint::measures
