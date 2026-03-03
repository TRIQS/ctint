// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3pp_iw_full.hpp"

namespace triqs_ctint::measures {

  M3pp_iw_full::M3pp_iw_full(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()), G0_tau(std::move(G0_tau_)) {

    // Construct full fermionic Matsubara mesh
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M3};

    // Init measurement container and capture view
    results->M3pp_iw_nfft_full = make_block2_gf(mesh::prod{iw_mesh, iw_mesh}, params.gf_struct);
    M3pp_iw_.rebind(results->M3pp_iw_nfft_full.value());
    M3pp_iw_() = 0;

    // Initialize intermediate scattering matrix on the same imfreq mesh (for type1 NFFT)
    GM = block_gf{iw_mesh, params.gf_struct};

    // Create type1 nfft buffers that write to the block_gf data
    for (auto bl : range(params.n_blocks())) {
      buf_arrarr(bl) = array_adapter{GM[bl].target_shape(), [&](int i, int j) {
        return nfft_buf_t{slice_target_to_scalar(GM[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
      }};
    }
  }

  void M3pp_iw_full::accumulate(mc_weight_t sign) {
    Z += sign;

    // Reset intermediate scattering matrix
    GM() = 0;

    // Init intermediate scattering matrices
    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = GM[bl].target_shape()[0];
      long k      = det.size();

      auto arr_GM = nda::zeros<dcomplex>(bl_size, bl_size, k);

      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (long j = 0; j < k; ++j) {
          auto &[tau_j, u_j, _, _, _] = det.get_y(j);
          auto G0_bj               = G0_tau[bl][closest_mesh_pt(params.beta - tau_j)];
          auto Ginv_ji             = det.inverse_matrix(j, i);
          for (int b_u : range(bl_size)) arr_GM(b_u, u_i, i) += G0_bj(b_u, u_j) * Ginv_ji;
        }
      }

      for (int m : range(bl_size))
        for (int n : range(bl_size))
          for (long i = 0; i < k; ++i) buf_arrarr(bl)(m, n).push_back({params.beta - det.get_x(i).tau}, arr_GM(m, n, i));

      for (auto &buf : buf_arrarr(bl)) buf.flush();
    }

    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {
        int bl1_size     = GM[bl1].target_shape()[0];
        int bl2_size     = GM[bl2].target_shape()[0];
        auto const &GM1  = GM[bl1];
        auto const &GM2  = GM[bl2];
        auto &M3pp_iw    = M3pp_iw_(bl1, bl2);

        // Loop over full frequency grid
        for (auto mp : M3pp_iw.mesh()) {
          auto [mp1, mp2] = mp;
          auto iw1 = mp1.value(); // matsubara_freq
          auto iw2 = mp2.value(); // matsubara_freq
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

  void M3pp_iw_full::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3pp_iw_ = mpi::all_reduce(M3pp_iw_, comm);
    M3pp_iw_ = M3pp_iw_ / Z;
  }

} // namespace triqs_ctint::measures
