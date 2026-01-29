// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3pp_iw.hpp"

namespace triqs_ctint::measures {

  M3pp_iw::M3pp_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()), G0_tau(std::move(G0_tau_)) {

    // Construct DLR2D Matsubara mesh
    mesh::dlr2d_imfreq M3pp_iw_mesh{params.beta, params.dlr_wmax_M3, params.dlr_eps_M3, mesh::PP};

    // Init measurement container and capture view
    results->M3pp_iw_nfft = make_block2_gf(M3pp_iw_mesh, params.gf_struct);
    M3pp_iw_.rebind(results->M3pp_iw_nfft.value());
    M3pp_iw_() = 0;

    // Initialize intermediate scattering matrix on regular imfreq mesh (for type1 NFFT)
    mesh::imfreq iw_mesh{params.beta, Fermion, M3pp_iw_mesh.max_n() + 1};
    GM = block_gf{iw_mesh, params.gf_struct};

    // Create type1 nfft buffers that write to the block_gf data
    for (auto bl : range(params.n_blocks())) {
      buf_arrarr(bl) = array_adapter{GM[bl].target_shape(), [&](int i, int j) {
        return nfft_buf_t{slice_target_to_scalar(GM[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
      }};
    }
  }

  void M3pp_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Reset intermediate scattering matrix
    GM() = 0;

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

    // Flush all buffers
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush();

    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        int bl1_size     = GM[bl1].target_shape()[0];
        int bl2_size     = GM[bl2].target_shape()[0];
        auto const &GM1  = GM[bl1];
        auto const &GM2  = GM[bl2];
        auto &M3pp_iw    = M3pp_iw_(bl1, bl2);

        // Single loop over DLR2D mesh points
        for (auto mp : M3pp_iw.mesh()) {
          auto [iw1, iw2] = mp.value(); // matsubara_freq pair
          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size)) {
                  M3pp_iw[mp](i, j, k, l) += sign * GM1[iw1](j, i) * GM2[iw2](l, k);
                  if (bl1 == bl2) { M3pp_iw[mp](i, j, k, l) -= sign * GM1[iw1](l, i) * GM2[iw2](j, k); }
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
