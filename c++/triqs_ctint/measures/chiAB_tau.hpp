// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "operator_block.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $\chi_{AB}(\tau)$ by operator insertion
  */
  struct chiAB_tau {

    chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    /// Accumulate chiAB_tau by operator insertion
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    params_t const &params;
    qmc_config_t &qmc_config;

    // Container for the accumulation
    gf_view<mesh::dlr_imtime, tensor_valued<1>> chiAB_tau_;

    // Unified operator term groups
    int n_blocks_;
    std::vector<term_group_t> groups_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Precomputed tau points from the DLR imtime mesh
    long L_;
    std::vector<tau_t> tau_points_;

    // Fixed tau for B-side operators (time 0)
    tau_t tau_zero_ = tau_t::zero();
  };

} // namespace triqs_ctint::measures
