// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  // Measure the histogram of perturbation order
  struct histogram {

    histogram(params_t const &, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate perturbation order into histogram
    void accumulate(mc_weight_t);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Reference to accumulation vector
    std::optional<std::vector<double>> &histogram_;

    // Accumulation counter
    long N = 0;
  };

} // namespace triqs_ctint::measures
