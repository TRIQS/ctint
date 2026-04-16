// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"
#include <triqs/stat/lin_binning.hpp>

namespace triqs_ctint::measures {

  /// Measure of the average perturbation order
  struct average_k {

    average_k(params_t const &, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate average perturbation order
    void accumulate(mc_weight_t);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    /// Report current value representation
    std::string report() const;

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Reference to double for accumulation
    double &average_k_;
    std::optional<double> &average_k_error_;

    // Linear binning for error estimation
    triqs::stat::lin_binning<dcomplex> k_bins_;

    // Accumulation counter
    long long N = 0;
  };

} // namespace triqs_ctint::measures
