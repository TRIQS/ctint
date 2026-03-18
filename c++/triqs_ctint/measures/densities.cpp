// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./densities.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  densities::densities(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), results_(results) {

    if (params.measure_densities) {
      results->densities = block_vector_t{};

      // Init measurement container and capture view, and init binning accumulators
      for (auto &[bl, bl_size] : params.gf_struct) {
        results->densities->push_back(nda::zeros<g_tau_scalar_t>(bl_size));
        densities_.push_back(results->densities->back());
        dens_bins_.emplace_back(nda::zeros<dcomplex>(bl_size), 128, 1);
      }
    }

    // Log-binning for auto-correlation: [0] = perturbation order, [1..] = sign * density per orbital
    int n_orbitals = 0;
    for (auto const &[bl, bl_size] : params.gf_struct) n_orbitals += bl_size;
    log_accs_.reserve(1 + n_orbitals);
    for (int i = 0; i < 1 + n_orbitals; ++i) log_accs_.emplace_back(dcomplex{0.0}, -1);
  }

  void densities::accumulate(mc_weight_t sign) {
    log_accs_[0] << dcomplex(qmc_config.perturbation_order());

    if (params.measure_densities) {
      Z += sign;
      ++N_;
    }

    // Measure diagonal densities using batched insert_ratios
    int idx = 1;
    for (int bl : range(params.n_blocks())) {

      auto &det   = qmc_config.dets[bl];
      int bl_size = params.gf_struct[bl].second;

      // Build paired c / cdag arrays for diagonal entries
      nda::array<c_t, 1> cs(bl_size);
      nda::array<cdag_t, 1> cdags(bl_size);
      for (int a = 0; a < bl_size; ++a) {
        cs(a)    = c_t{tau_t::get_zero(), a};
        cdags(a) = cdag_t{tau_t::get_zero_plus(), a};
      }

      // Single batched call: insert_ratios returns ratios for paired (cs[a], cdags[a])
      auto ratios = det.insert_ratios(0, 0, cs, cdags);
      for (int a = 0; a < bl_size; ++a) {
        auto val = sign * ratios(a);
        log_accs_[idx++] << val;
      }

      if (params.measure_densities) {
        auto step = nda::array<dcomplex, 1>(bl_size);
        for (int a = 0; a < bl_size; ++a) {
          auto val = sign * ratios(a);
          densities_[bl](a) += val;
          step(a) = val;
        }
        dens_bins_[bl] << step;
      }
    }
  }

  void densities::collect_results(mpi::communicator const &comm) {
    using triqs::stat::log_binning;

    // Auto-correlation time from log-binning
    results_->auto_corr_time = 0.0;
    for (auto &log_acc : log_accs_) {
      auto [mean, errs, taus, effs] = log_acc.mean_errors_and_taus(comm);
      if (!taus.empty()) { results_->auto_corr_time = std::max(results_->auto_corr_time, std::real(taus.back())); }
      log_acc = log_binning<dcomplex>{dcomplex{0.0}, -1};
    }
    mpi::broadcast(results_->auto_corr_time, comm, 0);

    // Densities
    if (params.measure_densities) {
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
  }

} // namespace triqs_ctint::measures
