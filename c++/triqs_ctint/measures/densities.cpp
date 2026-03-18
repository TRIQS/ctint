// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./densities.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  densities::densities(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), results_(results) {

    results->densities = block_vector_t{};

    // Init measurement container and capture view, and init binning accumulators
    for (auto &[bl, bl_size] : params.gf_struct) {
      results->densities->push_back(nda::zeros<g_tau_scalar_t>(bl_size));
      densities_.push_back(results->densities->back());
      dens_bins_.emplace_back(nda::zeros<dcomplex>(bl_size), 128, 1);
    }
  }

  void densities::accumulate(mc_weight_t sign) {
    Z += sign;
    ++N_;

    // Measure diagonal densities using batched insert_ratios
    for (int bl : range(params.n_blocks())) {

      auto &det   = qmc_config.dets[bl];
      auto &dens  = densities_[bl];
      int bl_size = dens.shape()[0];

      // Build paired c / cdag arrays for diagonal entries
      nda::array<c_t, 1> cs(bl_size);
      nda::array<cdag_t, 1> cdags(bl_size);
      for (int a = 0; a < bl_size; ++a) {
        cs(a)    = c_t{tau_t::get_zero(), a};
        cdags(a) = cdag_t{tau_t::get_zero_plus(), a};
      }

      // Single batched call: insert_ratios returns ratios for paired (cs[a], cdags[a])
      auto ratios = det.insert_ratios(0, 0, cs, cdags);
      auto step   = nda::array<dcomplex, 1>(bl_size);
      for (int a = 0; a < bl_size; ++a) {
        auto val = sign * ratios(a);
        dens(a) += val;
        step(a) = val;
      }
      dens_bins_[bl] << step;
    }
  }

  void densities::collect_results(mpi::communicator const &comm) {
    Z  = mpi::all_reduce(Z, comm);
    N_ = mpi::all_reduce(N_, comm);
    for (auto &dens : densities_) {
      dens = mpi::all_reduce(dens, comm);
      dens = dens / Z;
    }

    // Compute error bars from linear binning
    results_->densities_errors = block_vector_t{};
    for (int bl : range(params.n_blocks())) {
      auto [m, err, tau] = dens_bins_[bl].mean_error_and_tau(comm);
      auto norm          = std::abs(Z / N_);
      results_->densities_errors->push_back(nda::array<g_tau_scalar_t, 1>(nda::abs(err) / norm));
    }
  }

} // namespace triqs_ctint::measures
