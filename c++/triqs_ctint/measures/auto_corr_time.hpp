// Copyright (c) 2021--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/stat/accumulator.hpp>
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /// Measurement auto-correlation time based on the partition function
  struct auto_corr_time {

    auto_corr_time(params_t const &, qmc_config_t const &, container_set *results);

    /// Accumulate average sign
    void accumulate(mc_weight_t sign);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Reference to double for accumulation
    double &auto_corr_time_;

    // Initialize one complex log accumulator for each observable to use for the autocorrelation analysis
    std::vector<triqs::stat::accumulator<dcomplex>> log_accs = {2, {0.0, -1, 0}};
  };

} // namespace triqs_ctint::measures
