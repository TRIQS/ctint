// Copyright (c) 2023--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M4ph_iw.hpp"
#include <cmath>

namespace triqs_ctint::measures {

  M4ph_iw::M4ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    mesh::imfreq iW_mesh{params.beta, Boson, params.n_iW_M4};
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M4};
    mesh::prod<imfreq, imfreq, imfreq> M4ph_iw_mesh{iW_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M4ph_iw = make_block2_gf(M4ph_iw_mesh, params.gf_struct);
    M4ph_iw_.rebind(results->M4ph_iw.value());
    M4ph_iw_() = 0;

    // Construct Matsubara mesh for temporary Matrix
    mesh::imfreq iw_mesh_large{params.beta, Fermion, params.n_iW_M4 + params.n_iw_M4 + 1};
    mesh::prod<imfreq, imfreq> M_mesh{iw_mesh_large, iw_mesh_large};

    // Initialize intermediate scattering matrix
    M = block_gf{M_mesh, params.gf_struct};

    // Create nfft buffers
    for (int bl : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft_buf_t<2>{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta};
      };
      buf_arrarr(bl) = array_adapter{M[bl].target_shape(), init_target_func};
    }
  }

  void M4ph_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    M() = 0;
    for (int bl : range(params.n_blocks()))
      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[bl],
               [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) { // Care for negative frequency in c transform (for M-objects)
                 buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({double(cdag_j.tau), params.beta - double(c_i.tau)}, -Ginv_ji);
               });
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    auto const &iW_mesh = std::get<0>(M4ph_iw_(0, 0).mesh());
    auto const &iw_mesh = std::get<1>(M4ph_iw_(0, 0).mesh());

    for (int bl1 : range(params.n_blocks())) // FIXME c++17 Loops
      for (int bl2 : range(params.n_blocks())) {

        int bl1_size   = M[bl1].target_shape()[0];
        int bl2_size   = M[bl2].target_shape()[0];
        auto const &M1 = M[bl1];
        auto const &M2 = M[bl2];
        auto &M4       = M4ph_iw_(bl1, bl2);

        for (auto iW : iW_mesh)
          for (auto iw : iw_mesh)
            for (auto iwp : iw_mesh)
              for (int i : range(bl1_size))
                for (int j : range(bl1_size)) {
                  auto M1val = M1[iW + iw, iw.value()](j, i);

                  for (int k : range(bl2_size))
                    for (int l : range(bl2_size)) {
                      M4[iW, iw, iwp](i, j, k, l) += sign * M1val * M2[iwp.value(), iW + iwp](l, k);
                      if (bl1 == bl2) { M4[iW, iw, iwp](i, j, k, l) -= sign * M1[iwp.value(), iw.value()](l, i) * M2[iW + iw, iW + iwp](j, k); }
                    }
                }
      }
  }

  void M4ph_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M4ph_iw_ = mpi::all_reduce(M4ph_iw_, comm);
    M4ph_iw_ = M4ph_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
