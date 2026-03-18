// Copyright (c) 2021--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/stat/log_binning.hpp>
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /// Measurement of the auto-correlation time based on perturbation order and diagonal densities
  struct auto_corr_time {

    auto_corr_time(params_t const &params, qmc_config_t const &qmc_config, container_set *results);

    /// Accumulate observables for autocorrelation analysis
    void accumulate(mc_weight_t sign);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    params_t const &params;
    qmc_config_t const &qmc_config;
    double &auto_corr_time_;

    // Log-binning accumulators: [0] = perturbation order, [1..] = sign * density per orbital
    std::vector<triqs::stat::log_binning<dcomplex>> log_accs;
  };

} // namespace triqs_ctint::measures
