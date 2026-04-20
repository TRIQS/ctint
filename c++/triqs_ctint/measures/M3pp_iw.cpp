// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M3pp_iw.hpp"
#include "./iw_accumulate.hpp"

namespace triqs_ctint::measures {

  M3pp_iw::M3pp_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()), G0_tau(std::move(G0_tau_)) {

    // Construct DLR2D Matsubara mesh
    mesh::dlr2d_imfreq M3pp_iw_mesh{params.beta, params.dlr_wmax, params.dlr_eps, mesh::PP, params.dlr2d_compress_grid};

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
        return nfft::buffer_t{slice_target_to_scalar(GM[bl], i, j).data(), params.nfft_buf_size, params.beta, params.nfft_tol};
      }};
    }
  }

  void M3pp_iw::accumulate(mc_weight_t sign) {
    Z += sign;

    // Reset intermediate scattering matrix
    GM() = 0;

    // Compute GM = G0 * M^{-1} via BLAS GEMM instead of O(k^2 * bl_size) scalar loop
    for (int bl : range(params.n_blocks())) {
      auto &det   = qmc_config.dets[bl];
      int bl_size = GM[bl].target_shape()[0];
      long k      = det.size();
      if (k == 0) continue;

      // Build Y(b, j) = G0(beta - tau_j)(b, u_j) — one G0 lookup per j instead of per (i,j)
      auto Y = nda::matrix<dcomplex>(bl_size, k);
      for (long j = 0; j < k; ++j) {
        auto &[tau_j, u_j, _, _, _] = det.get_y(j);
        auto G0_j                   = G0_tau[bl][closest_mesh_pt(params.beta - tau_j)];
        for (int b = 0; b < bl_size; ++b) Y(b, j) = G0_j(b, u_j);
      }

      // W = Y * Ginv via BLAS GEMM: (bl_size, k) * (k, k) -> (bl_size, k)
      auto W = Y * det.inverse_matrix();

      // Push non-zero entries to NFFT buffers (skip orbital indices with zero contribution)
      for (long i = 0; i < k; ++i) {
        auto &[tau_i, u_i, _, _, _] = det.get_x(i);
        for (int m : range(bl_size)) buf_arrarr(bl)(m, u_i).push_back({params.beta - tau_i}, W(m, i));
      }

      for (auto &buf : buf_arrarr(bl)) buf.flush();
    }

    for (auto bl1 : range(params.n_blocks()))
      for (auto bl2 : range(params.n_blocks())) {
        auto const bl2_size = GM[bl2].target_shape()[0];
        simd::dlr2d_iw3pp_accumulate(sign, GM, M3pp_iw_, bl1, bl2, bl2_size);
      }
  }

  void M3pp_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M3pp_iw_ = mpi::all_reduce(M3pp_iw_, comm);
    M3pp_iw_ = M3pp_iw_ / Z;
  }

} // namespace triqs_ctint::measures
