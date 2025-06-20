// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./density.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  density::density(params_t const &params_, qmc_config_t &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    results->density = block_matrix_t{};

    // Init measurement container and capture view
    for (auto &[bl, bl_size] : params.gf_struct) {
      results->density->push_back(zeros<M_tau_scalar_t>(make_shape(bl_size, bl_size)));
      density_.push_back(results->density->back());
    }
  }

  void density::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate density matrix
    for (int bl : range(params.n_blocks())) {

      auto &det     = qmc_config.dets[bl];
      auto &density = density_[bl];
      int bl_size   = density.shape()[0];

      for (int a : range(bl_size))
        for (int b : range(bl_size)) {

          cdag_t cdag_a{tau_t::get_zero_plus(), a};
          c_t c_b{tau_t::get_zero(), b};

          density(a, b) += sign * det.try_insert(0, 0, c_b, cdag_a);
          det.reject_last_try();
        }
    }
  }

  void density::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z = mpi::all_reduce(Z, comm);
    for (auto &dens_mat : density_) {
      dens_mat = mpi::all_reduce(dens_mat, comm);
      dens_mat = dens_mat / Z;
    }
  }

} // namespace triqs_ctint::measures
