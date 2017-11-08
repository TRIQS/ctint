// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /// Measure of the average sign
  struct average_sign {

    average_sign(params_t const &, qmc_config_t const &, container_set *results);

    /// Accumulate average sign
    void accumulate(mc_weight_t sign);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    /// Report the current value representation
    std::string report() const;

    private:
    // Reference to double for accumulation
    mc_weight_t &average_sign_, average_static_sign, average_dynamic_sign;
    uint64_t &nmeasures;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Accumulation counter
    long count = 0;
  };

} // namespace triqs_ctint::measures
