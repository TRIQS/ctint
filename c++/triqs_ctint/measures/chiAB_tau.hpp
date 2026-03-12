// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include "../nfft/buffer.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $\chi_{AB}(\tau)$ by operator insertion
  */
  struct chiAB_tau {

    chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters FIXME We cannot choose const, as we call try_insert
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t &qmc_config;

    // Container for the accumulation
    gf_view<mesh::dlr_imtime, tensor_valued<1>> chiAB_tau_;

    // The parsed operator pairs
    using op_term_t = std::tuple<dcomplex, std::pair<int, int>, std::pair<int, int>>;
    std::vector<std::pair<std::vector<op_term_t>, std::vector<op_term_t>>> op_pairs;

    // The average sign
    mc_weight_t Z = 0.0;

    // The tau-mesh
    mesh::dlr_imtime tau_mesh;
  };

} // namespace triqs_ctint::measures
