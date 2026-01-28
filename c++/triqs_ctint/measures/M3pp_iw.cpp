// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3pp_iw.hpp"
#include "./M3_iw_utils.hpp"
#include <set>

namespace triqs_ctint::measures {

  M3pp_iw::M3pp_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()), G0_tau(std::move(G0_tau_)) {

    // Construct DLR2D Matsubara mesh
    mesh::dlr2d_imfreq M3pp_iw_mesh{params.beta, params.dlr_wmax_M3, params.dlr_eps_M3, mesh::PP};

    // Init measurement container and capture view
    results->M3pp_iw_nfft = make_block2_gf(M3pp_iw_mesh, params.gf_struct);
    M3pp_iw_.rebind(results->M3pp_iw_nfft.value());
    M3pp_iw_() = 0;

    // Collect all unique Matsubara frequencies from DLR2D mesh
    std::set<mesh::matsubara_freq> unique_w_set;
    for (auto [w1, w2] : M3pp_iw_mesh) {
      unique_w_set.insert(w1);
      unique_w_set.insert(w2);
    }
    // Build target matsubara_freq vector and index map
    std::tie(target_mf_1d, n_idx_offset, n_to_idx) = build_index_map(unique_w_set);

    // Initialize GM_data arrays and create type 3 NFFT buffers
    GM_data.resize(params.n_blocks());
    for (auto bl : range(params.n_blocks())) {
      auto bl_size = params.gf_struct[bl].second;
      GM_data(bl).resize(target_mf_1d.size(), bl_size, bl_size);
      GM_data(bl) = 0;

      buf_arrarr(bl) = array_adapter{std::array{bl_size, bl_size}, [&](int i, int j) {
        return nfft_buf_t{GM_data(bl)(nda::range::all, i, j), target_mf_1d, params.nfft_buf_size, nfft_type_t::type3, params.nfft_tol};
      }};
    }
  }

  void M3pp_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Reset intermediate scattering matrix
    for (auto &gm : GM_data) gm = 0;

    // Fill GM using explicit row-column loops to accumulate values before pushing to NFFT buffers:
    // For fixed c_i (row), sum over all cdag_j (columns) before pushing once
    for (int bl : range(params.n_blocks())) {
      int bl_size = params.gf_struct[bl].second;
      auto &det   = qmc_config.dets[bl];
      long k_bl   = det.size();

      nda::array<dcomplex, 1> gm_acc(bl_size);

      for (long row = 0; row < k_bl; ++row) {
        auto const &c_i = det.get_x(row);
        auto tau_i      = double(c_i.tau);

        gm_acc = 0;
        for (long col = 0; col < k_bl; ++col) {
          auto const &cdag_j = det.get_y(col);
          auto Ginv_ji       = det.inverse_matrix(col, row);
          auto G0_btau_j     = G0_tau[bl][closest_mesh_pt(params.beta - double(cdag_j.tau))];
          for (int b_u : range(bl_size)) gm_acc(b_u) += G0_btau_j(b_u, cdag_j.u) * Ginv_ji;
        }

        // Push accumulated GM once per (b_u, c_i.u)
        for (int b_u : range(bl_size)) buf_arrarr(bl)(b_u, c_i.u).push_back({params.beta - tau_i}, gm_acc(b_u));
      }
    }
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush();

    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        int bl1_size  = params.gf_struct[bl1].second;
        int bl2_size  = params.gf_struct[bl2].second;
        auto &M3pp_iw = M3pp_iw_(bl1, bl2);

        // Single loop over DLR2D mesh points
        for (auto mp : M3pp_iw.mesh()) {
          auto [n1, n2] = mp.index();
          auto idx1     = n_to_idx[n1 + n_idx_offset];
          auto idx2     = n_to_idx[n2 + n_idx_offset];
          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size)) {
                  M3pp_iw[mp](i, j, k, l) += sign * GM_data(bl1)(idx1, j, i) * GM_data(bl2)(idx2, l, k);
                  if (bl1 == bl2) { M3pp_iw[mp](i, j, k, l) -= sign * GM_data(bl1)(idx1, l, i) * GM_data(bl2)(idx2, j, k); }
                }
        }
      }
  }

  void M3pp_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3pp_iw_ = mpi::all_reduce(M3pp_iw_, comm);
    M3pp_iw_ = M3pp_iw_ / Z;
  }

} // namespace triqs_ctint::measures
