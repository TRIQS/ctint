// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3ph_iw_full.hpp"

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
                         return nfft_buf_t{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
                       }};
      buf_arrarr_GM(bl) =
         array_adapter{GM[bl].target_shape(), [&](int i, int j) {
                         return nfft_buf_t{slice_target_to_scalar(GM[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
                       }};
      buf_arrarr_MG(bl) =
         array_adapter{MG[bl].target_shape(), [&](int i, int j) {
                         return nfft_buf_t{slice_target_to_scalar(MG[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
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

    // Init intermediate scattering matrices
    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = GM[bl].target_shape()[0];
      long k      = det.size();

      auto arr_GM = nda::zeros<dcomplex>(bl_size, bl_size, k);
      auto arr_MG = nda::zeros<dcomplex>(bl_size, bl_size, k);

      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (long j = 0; j < k; ++j) {
          auto &[tau_j, u_j, _, _, _] = det.get_y(j);

          auto Ginv_ji = det.inverse_matrix(j, i);

          // Fill M, Note: Minus sign from the shift of -tau_i
          buf_arrarr(bl)(u_j, u_i).push_back({tau_j, beta - tau_i}, -Ginv_ji);

          //Fill GMG, GM, MG
          for (int abar_u : range(bl_size)) {
            auto G0_ia = G0_tau[bl][closest_mesh_pt(double(tau_i))](u_i, abar_u);
            for (int b_u : range(bl_size)) {
              // Note: Minus sign from the shift of -tau_j
              auto G0_bj = -G0_tau[bl][closest_mesh_pt(beta - tau_j)](b_u, u_j);
              GMG(bl)(abar_u, b_u) += G0_bj * Ginv_ji * G0_ia;
              // Note: Minus sign from the shift of -tau_i
              arr_GM(b_u, u_i, i) += -G0_bj * Ginv_ji;
              arr_MG(b_u, u_i, j) += Ginv_ji * G0_ia;
            }
          }
        }
      }
      for (auto m : range(bl_size)) {
        for (auto n : range(bl_size)) {
          for (long p = 0; p < k; ++p) {
            buf_arrarr_GM(bl)(m, n).push_back({beta - det.get_x(p).tau}, arr_GM(m, n, p));
            buf_arrarr_MG(bl)(m, n).push_back({det.get_y(p).tau}, arr_MG(m, n, p));
          }
        }
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

        // Loop over full frequency grid
        for (auto mp : M3ph_iw.mesh()) {
          auto [mp1, mp2] = mp;
          auto iw1 = mp1.value(); // matsubara_freq
          auto iw2 = mp2.value(); // matsubara_freq
          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size)) {
                  // M1 is on same prod mesh; access transposed (iw2, iw1)
                  M3ph_iw[mp](i, j, k, l) += sign * M1[closest_mesh_pt(iw2, iw1)](j, i) * GMG2(l, k);
                  if (bl1 == bl2) { M3ph_iw[mp](i, j, k, l) -= sign * GM1[iw1](l, i) * MG2[iw2](j, k); }
                }
        }
      }
  }

  void M3ph_iw_full::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3ph_iw_ = mpi::all_reduce(M3ph_iw_, comm);
    M3ph_iw_ = M3ph_iw_ / Z;
  }

} // namespace triqs_ctint::measures
