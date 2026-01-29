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
    mesh::dlr2d_imfreq M3ph_iw_mesh{params.beta, params.dlr_wmax, params.dlr_eps, mesh::PH};

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
    Z += sign;

    // Reset all accumulators
    for (auto &gmg : GMG) gmg() = 0;
    M() = 0;
    GM() = 0;
    MG() = 0;

    double beta = params.beta;

    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = params.gf_struct[bl].second;
      long k      = det.size();

      // ═══════════════════════════════════════════════════════════════════════
      // Phase 1: M(j.u, i.u)[iw1, iw2] — 2D NFFT at (tau_j, beta-tau_i)
      // Cannot reduce push_backs: each (i,j) pair has unique (tau_i, tau_j)
      // ═══════════════════════════════════════════════════════════════════════
      for (long row = 0; row < k; ++row) {
        auto const &c_i = det.get_x(row);
        double tau_i    = double(c_i.tau);
        for (long col = 0; col < k; ++col) {
          auto const &cdag_j = det.get_y(col);
          // Note: Minus sign from the shift of -tau_i
          buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({double(cdag_j.tau), beta - tau_i}, -det.inverse_matrix(col, row));
        }
      }

      // ═══════════════════════════════════════════════════════════════════════
      // Phase 2: GMG(a, b) — no NFFT, direct accumulation
      // GMG(a,b) = sum_{i,j} G0(beta-tau_j)(b,j.u) * Ginv(j,i) * G0(tau_i)(i.u,a)
      // ═══════════════════════════════════════════════════════════════════════
      for (long row = 0; row < k; ++row) {
        auto const &c_i  = det.get_x(row);
        auto G0_at_tau_i = G0_tau[bl][closest_mesh_pt(double(c_i.tau))];
        for (long col = 0; col < k; ++col) {
          auto const &cdag_j   = det.get_y(col);
          auto G0_at_btau_j    = G0_tau[bl][closest_mesh_pt(beta - double(cdag_j.tau))];
          auto Ginv_ji         = det.inverse_matrix(col, row);
          for (int a : range(bl_size))
            for (int b : range(bl_size))
              // Note: Minus sign from the shift of -tau_j
              GMG(bl)(a, b) += (-G0_at_btau_j(b, cdag_j.u)) * Ginv_ji * G0_at_tau_i(c_i.u, a);
        }
      }

      // ═══════════════════════════════════════════════════════════════════════
      // Phase 3: GM(b, i.u)[iw] — 1D NFFT at (beta-tau_i)
      // For fixed row i: GM(b, i.u) = sum_j G0(beta-tau_j)(b, j.u) * Ginv(j,i)
      // Original code had factor bl_size from sum over abar_u index
      // ═══════════════════════════════════════════════════════════════════════
      nda::array<dcomplex, 1> gm_acc(bl_size);
      for (long row = 0; row < k; ++row) {
        auto const &c_i = det.get_x(row);
        gm_acc          = 0;
        for (long col = 0; col < k; ++col) {
          auto const &cdag_j = det.get_y(col);
          auto G0_at_btau_j  = G0_tau[bl][closest_mesh_pt(beta - double(cdag_j.tau))];
          auto Ginv_ji       = det.inverse_matrix(col, row);
          for (int b : range(bl_size))
            // Note: Minus sign from the shift of -tau_j
            gm_acc(b) += (-G0_at_btau_j(b, cdag_j.u)) * Ginv_ji;
        }
        // Push once per row; bl_size factor from original abar_u summation
        // Note: Minus sign from the shift of -tau_i
        for (int b : range(bl_size))
          buf_arrarr_GM(bl)(b, c_i.u).push_back({beta - double(c_i.tau)}, -double(bl_size) * gm_acc(b));
      }

      // ═══════════════════════════════════════════════════════════════════════
      // Phase 4: MG(b, i.u)[iw] — 1D NFFT at (tau_j)
      // For fixed col j: MG(b, i.u) = sum_i Ginv(j,i) * sum_a G0(tau_i)(i.u, a)
      // ═══════════════════════════════════════════════════════════════════════
      nda::array<dcomplex, 2> mg_acc(k, bl_size);
      mg_acc = 0;
      for (long row = 0; row < k; ++row) {
        auto const &c_i  = det.get_x(row);
        auto G0_at_tau_i = G0_tau[bl][closest_mesh_pt(double(c_i.tau))];
        dcomplex g0_row_sum = 0;
        for (int a : range(bl_size)) g0_row_sum += G0_at_tau_i(c_i.u, a);
        for (long col = 0; col < k; ++col)
          mg_acc(col, c_i.u) += det.inverse_matrix(col, row) * g0_row_sum;
      }
      for (long col = 0; col < k; ++col) {
        double tau_j = double(det.get_y(col).tau);
        for (int u_i : range(bl_size))
          for (int b : range(bl_size))
            buf_arrarr_MG(bl)(b, u_i).push_back({tau_j}, mg_acc(col, u_i));
      }
    }

    // Flush all buffers
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
