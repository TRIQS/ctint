// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3ph_iw.hpp"
#include "./M3_iw_utils.hpp"
#include <set>

namespace triqs_ctint::measures {

  M3ph_iw::M3ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_),
       qmc_config(qmc_config_),
       buf_arrarr(params_.n_blocks()),
       buf_arrarr_GM(params_.n_blocks()),
       buf_arrarr_MG(params_.n_blocks()),
       G0_tau(std::move(G0_tau_)),
       M_data(params.n_blocks()),
       GM_data(params.n_blocks()),
       MG_data(params.n_blocks()),
       GMG(params.n_blocks()) {

    // Construct DLR2D Matsubara mesh
    mesh::dlr2d_imfreq M3ph_iw_mesh{params.beta, params.dlr_wmax_M3, params.dlr_eps_M3, mesh::PH};
    int64_t n_mesh_points = M3ph_iw_mesh.size();

    // Init measurement container and capture view
    results->M3ph_iw_nfft = make_block2_gf(M3ph_iw_mesh, params.gf_struct);
    M3ph_iw_.rebind(results->M3ph_iw_nfft.value());
    M3ph_iw_() = 0;

    // Collect unique w1 (for GM) and w2 (for MG) and build 2D target for M
    std::set<mesh::matsubara_freq> unique_w1_set, unique_w2_set;
    // push_back uses {tau_j, beta-tau_i}: dim 0 targets w2, dim 1 targets w1
    target_mf_2d.reserve(n_mesh_points);
    for (auto [w1, w2] : M3ph_iw_mesh) {
      unique_w1_set.insert(w1);
      unique_w2_set.insert(w2);
      target_mf_2d.push_back({w2, w1});
    }
    // Build 1D target matsubara_freq vectors and index maps for GM (w1) and MG (w2)
    std::tie(target_mf_n1, n1_idx_offset, n1_to_idx) = build_index_map(unique_w1_set);
    std::tie(target_mf_n2, n2_idx_offset, n2_to_idx) = build_index_map(unique_w2_set);

    // Initialize data arrays, GMG, and NFFT buffers

    for (auto bl : range(params.n_blocks())) {
      auto bl_size = params.gf_struct[bl].second;
      GMG(bl)      = array<dcomplex, 2>(bl_size, bl_size);

      M_data(bl).resize(n_mesh_points, bl_size, bl_size);
      GM_data(bl).resize(target_mf_n1.size(), bl_size, bl_size);
      MG_data(bl).resize(target_mf_n2.size(), bl_size, bl_size);

      buf_arrarr(bl) = array_adapter{std::array{bl_size, bl_size}, [&](int i, int j) {
        return nfft_buf_t{M_data(bl)(nda::range::all, i, j), target_mf_2d, params.nfft_buf_size, nfft_type_t::type3, params.nfft_tol};
      }};
      buf_arrarr_GM(bl) = array_adapter{std::array{bl_size, bl_size}, [&](int i, int j) {
        return nfft_buf_t{GM_data(bl)(nda::range::all, i, j), target_mf_n1, params.nfft_buf_size, nfft_type_t::type3, params.nfft_tol};
      }};
      buf_arrarr_MG(bl) = array_adapter{std::array{bl_size, bl_size}, [&](int i, int j) {
        return nfft_buf_t{MG_data(bl)(nda::range::all, i, j), target_mf_n2, params.nfft_buf_size, nfft_type_t::type3, params.nfft_tol};
      }};
    }
  }

  void M3ph_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Reset intermediate scattering matrices
    for (auto &gmg : GMG) { gmg() = 0; }
    for (auto &m : M_data) m = 0;
    for (auto &gm : GM_data) gm = 0;
    for (auto &mg : MG_data) mg = 0;

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
          auto G0_btau_j     = G0_tau[bl][closest_mesh_pt(beta - tau_j)];

          // Fill M, Note: Minus sign from the shift of -tau_i
          buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({tau_j, beta - tau_i}, -Ginv_ji);

          // Accumulate GM: sum_j G0(beta-tau_j)(b_u, u_j) * Ginv_ji
          for (int b_u : range(bl_size)) gm_acc(b_u) += G0_btau_j(b_u, cdag_j.u) * Ginv_ji;

          // Accumulate MG: sum over rows for this column and orbital
          mg_acc(col, c_i.u) += Ginv_ji * mg_g0_sum;

          // Fill GMG
          for (int abar_u : range(bl_size)) {
            auto G0_ia = G0_tau_i(c_i.u, abar_u);
            for (int b_u : range(bl_size))
              // Note: Minus sign from the shift of -tau_j
              GMG(bl)(abar_u, b_u) += -G0_btau_j(b_u, cdag_j.u) * Ginv_ji * G0_ia;
          }
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

        int bl1_size     = params.gf_struct[bl1].second;
        int bl2_size     = params.gf_struct[bl2].second;
        auto const &GMG2 = GMG(bl2);
        auto &M3ph_iw    = M3ph_iw_(bl1, bl2);

        // Single loop over DLR2D mesh points
        for (auto mp : M3ph_iw.mesh()) {
          auto [n1, n2] = mp.index();
          auto gm_idx   = n1_to_idx[n1 + n1_idx_offset];
          auto mg_idx   = n2_to_idx[n2 + n2_idx_offset];
          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size)) {
                  M3ph_iw[mp](i, j, k, l) += sign * M_data(bl1)(mp.data_index(), j, i) * GMG2(l, k);
                  if (bl1 == bl2) { M3ph_iw[mp](i, j, k, l) -= sign * GM_data(bl1)(gm_idx, l, i) * MG_data(bl2)(mg_idx, j, k); }
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
