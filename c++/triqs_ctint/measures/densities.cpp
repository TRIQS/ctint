// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./densities.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  densities::densities(params_t const &params_, qmc_config_t &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    results->densities = block_vector_t{};

    // Init measurement container and capture view
    for (auto &[bl, bl_size] : params.gf_struct) {
      results->densities->push_back(nda::zeros<g_tau_scalar_t>(bl_size));
      densities_.push_back(results->densities->back());
    }
  }

  void densities::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Measure diagonal densities using batched insert_ratios
    for (int bl : range(params.n_blocks())) {

      auto &det  = qmc_config.dets[bl];
      auto &dens = densities_[bl];
      int bl_size = dens.shape()[0];

      // Build paired c / cdag vectors for diagonal entries
      std::vector<c_t> cs(bl_size);
      std::vector<cdag_t> cdags(bl_size);
      for (int a = 0; a < bl_size; ++a) {
        cs[a]    = c_t{tau_t::get_zero(), a};
        cdags[a] = cdag_t{tau_t::get_zero_plus(), a};
      }

      // Single batched call: insert_ratios returns ratios for paired (cs[a], cdags[a])
      auto ratios = det.insert_ratios(0, 0, cs, cdags);
      for (int a = 0; a < bl_size; ++a) dens(a) += sign * ratios(a);
    }
  }

  void densities::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z = mpi::all_reduce(Z, comm);
    for (auto &dens : densities_) {
      dens = mpi::all_reduce(dens, comm);
      dens = dens / Z;
    }
  }

} // namespace triqs_ctint::measures
