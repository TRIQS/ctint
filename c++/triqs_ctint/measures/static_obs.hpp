// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/stat/lin_binning.hpp>
#include "operator_block.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of static expectation values $\langle C_i \rangle$ by operator insertion with tau-averaging.
  * Each C_i is a many_body_operator containing monomials of arbitrary even degree.
  * Time-translational invariance is exploited: C is evaluated at L uniformly spaced tau points and averaged.
  */
  struct static_obs {

    static_obs(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    void accumulate(mc_weight_t sign);

    void collect_results(mpi::communicator const &comm);

    private:
    params_t const &params;
    qmc_config_t &qmc_config;

    // Result storage
    container_set *results_;
    nda::array<dcomplex, 1> result_;

    // Unified operator term groups
    int n_blocks_;
    std::vector<term_group_t> groups_;

    // Constant contributions (degree-0 monomials)
    nda::array<dcomplex, 1> constant_parts_;

    // Tau averaging
    long L_;                        // number of tau points
    std::vector<tau_t> tau_points_; // precomputed uniform tau grid

    mc_weight_t Z = 0.0;
    long N_       = 0; // step counter

    // Per-step contributions and linear binning for error analysis
    nda::array<dcomplex, 1> step_contrib_;
    std::vector<triqs::stat::lin_binning<dcomplex>> obs_bins_;
  };

} // namespace triqs_ctint::measures
