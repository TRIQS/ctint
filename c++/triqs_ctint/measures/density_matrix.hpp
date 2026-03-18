// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/stat/lin_binning.hpp>
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of the full density matrix by operator insertion using insert_ratios_matrix
  */
  struct density_matrix {

    density_matrix(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    /// Accumulate density matrix using batched insert_ratios_matrix
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    params_t const &params;
    qmc_config_t &qmc_config;
    container_set *results_;

    // Container for the accumulation (views into results->density_matrix)
    block_matrix_v_t density_matrix_;

    mc_weight_t Z = 0.0;
    long N_        = 0;

    // Linear binning for error analysis (one accumulator per block)
    std::vector<triqs::stat::lin_binning<nda::array<dcomplex, 2>>> dm_bins_;
  };

} // namespace triqs_ctint::measures
