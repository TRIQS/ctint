// Copyright (c) 2021--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./auto_corr_time.hpp"

namespace triqs_ctint::measures {

  auto_corr_time::auto_corr_time(params_t const &params, qmc_config_t const &qmc_config, container_set *results)
     : params(params), qmc_config(qmc_config), auto_corr_time_(results->auto_corr_time) {

    // Count total density orbitals
    int n_orbitals = 0;
    for (auto const &[bl, bl_size] : params.gf_struct) n_orbitals += bl_size;

    // [0] = perturbation order, [1..] = sign * density per orbital
    log_accs.reserve(1 + n_orbitals);
    for (int i = 0; i < 1 + n_orbitals; ++i) log_accs.emplace_back(dcomplex{0.0}, -1);
  }

  void auto_corr_time::accumulate(mc_weight_t sign) {
    log_accs[0] << dcomplex(qmc_config.perturbation_order());

    // Accumulate sign * diagonal densities
    int idx = 1;
    for (int bl = 0; bl < params.n_blocks(); ++bl) {
      auto const &det = qmc_config.dets[bl];
      int bl_size     = params.gf_struct[bl].second;

      nda::array<c_t, 1> cs(bl_size);
      nda::array<cdag_t, 1> cdags(bl_size);
      for (int a = 0; a < bl_size; ++a) {
        cs(a)    = c_t{tau_t::get_zero(), a};
        cdags(a) = cdag_t{tau_t::get_zero_plus(), a};
      }

      auto ratios = det.insert_ratios(0, 0, cs, cdags);
      for (int a = 0; a < bl_size; ++a) log_accs[idx++] << sign * ratios(a);
    }
  }

  void auto_corr_time::collect_results(mpi::communicator const &comm) {
    using triqs::stat::log_binning;

    auto_corr_time_ = 0.0;

    for (auto &log_acc : log_accs) {
      auto [mean, errs, taus, effs] = log_acc.mean_errors_and_taus(comm);
      if (!taus.empty()) { auto_corr_time_ = std::max(auto_corr_time_, std::real(taus.back())); }

      // Reset the accumulator
      log_acc = log_binning<dcomplex>{dcomplex{0.0}, -1};
    }
    mpi::broadcast(auto_corr_time_, comm, 0);
  }

} // namespace triqs_ctint::measures
