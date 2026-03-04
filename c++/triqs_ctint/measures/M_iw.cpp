// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M_iw.hpp"

namespace triqs_ctint::measures {

  M_iw::M_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), M_data(params_.n_blocks()) {

    // Construct DLR Matsubara mesh
    mesh::dlr_imfreq M_iw_mesh{params.beta, Fermion, params.dlr_wmax, params.dlr_eps, true};
    int64_t n_dlr_pts = M_iw_mesh.size();

    // Init measurement container and capture view
    results->M_iw_nfft = g_dlr_iw_t{M_iw_mesh, params.gf_struct};
    M_iw_.rebind(results->M_iw_nfft.value());
    M_iw_() = 0;

    // Build target matsubara_freq vector
    target_mf.reserve(n_dlr_pts);
    for (auto w : M_iw_mesh) target_mf.push_back(w);

    // Initialize M_data arrays and create type 3 NFFT buffers
    for (int bl : range(params.n_blocks())) {
      int bl_size = params.gf_struct[bl].second;
      M_data(bl).resize(n_dlr_pts, bl_size, bl_size);
      M_data(bl) = 0;

      auto init_func = [&](int i, int j) {
        return nfft::buffer_t<1>{M_data(bl)(nda::range::all, i, j), target_mf, params.nfft_buf_size, params.nfft_tol};
      };
      buf_vec.emplace_back(array_adapter{std::array{bl_size, bl_size}, init_func});
    }

    // Initialize M_hartree if not already set (e.g. by M_tau measurement)
    if (!results->M_hartree) {
      results->M_hartree = make_block_vector<g_tau_scalar_t>(params.gf_struct);
      for (auto &m : results->M_hartree.value()) M_hartree_.push_back(m);
    }
  }

  void M_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Reset intermediate data
    for (auto &m : M_data) m = 0;

    // Loop over blocks
    for (int b = 0; b < M_iw_.size(); ++b) {
      // Loop over every index pair (x,y) in the determinant matrix
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {
        // Handle equal-time case (Hartree term) separately
        if (c_i.tau == cdag_j.tau) {
          if (!M_hartree_.empty()) M_hartree_[b](cdag_j.u, c_i.u) += Ginv * sign;
        } else {
          // Absolut time-difference tau of the index pair
          auto [s, dtau] = cyclic_difference(cdag_j.tau, c_i.tau);

          // Push {tau, f(tau)} pair into nfft buffer
          auto &buf = buf_vec[b](cdag_j.u, c_i.u);
          buf.push_back({dtau}, Ginv * s * sign);
        }
      });
    }

    // Flush buffers and copy to M_iw_
    for (auto &buf_arr : buf_vec)
      for (auto &buf : buf_arr) buf.flush();

    // Copy M_data to M_iw_
    for (int bl : range(params.n_blocks())) {
      int bl_size = params.gf_struct[bl].second;
      auto &M_bl  = M_iw_[bl];
      int64_t idx = 0;
      for (auto mp : M_bl.mesh()) {
        for (int i : range(bl_size))
          for (int j : range(bl_size)) M_bl[mp](i, j) += M_data(bl)(idx, i, j);
        ++idx;
      }
    }
  }

  void M_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z     = mpi::all_reduce(Z, comm);
    M_iw_ = mpi::all_reduce(M_iw_, comm);
    M_iw_ = M_iw_ / (-Z * params.beta);

    // Normalize M_hartree (only if M_iw is responsible for it)
    for (auto &m : M_hartree_) {
      m = mpi::all_reduce(m, comm);
      m = m / (-Z * params.beta);
    }
  }

} // namespace triqs_ctint::measures
