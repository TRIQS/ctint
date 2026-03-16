// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./density_matrix.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  density_matrix::density_matrix(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_) {

    results->density_matrix = block_matrix_t{};

    // Init measurement container and capture view
    for (auto &[bl, bl_size] : params.gf_struct) {
      results->density_matrix->push_back(zeros<g_tau_scalar_t>(make_shape(bl_size, bl_size)));
      density_matrix_.push_back(results->density_matrix->back());
    }
  }

  void density_matrix::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Measure full density matrix using batched insert_ratios_matrix
    for (int bl : range(params.n_blocks())) {

      auto &det      = qmc_config.dets[bl];
      auto &dens_mat = density_matrix_[bl];
      int bl_size    = dens_mat.shape()[0];

      // Build c and cdag vectors for all orbital indices
      std::vector<c_t> cs(bl_size);
      std::vector<cdag_t> cdags(bl_size);
      for (int a = 0; a < bl_size; ++a) {
        cs[a]    = c_t{tau_t::get_zero(), a};
        cdags[a] = cdag_t{tau_t::get_zero_plus(), a};
      }

      // Single batched call: insert_ratios_matrix returns ratios(b,a) for (c_b, cdag_a)
      // density_matrix(a,b) = <cdag_a c_b> = det_ratio for inserting (c_b, cdag_a)
      auto ratios = det.insert_ratios_matrix(0, 0, cs, cdags);
      // ratios(b, a) = det_ratio for (c_b, cdag_a) = density_matrix(a,b)
      for (int a = 0; a < bl_size; ++a)
        for (int b = 0; b < bl_size; ++b) dens_mat(a, b) += sign * ratios(b, a);
    }
  }

  void density_matrix::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z = mpi::all_reduce(Z, comm);
    for (auto &dens_mat : density_matrix_) {
      dens_mat = mpi::all_reduce(dens_mat, comm);
      dens_mat = dens_mat / Z;
    }
  }

} // namespace triqs_ctint::measures
