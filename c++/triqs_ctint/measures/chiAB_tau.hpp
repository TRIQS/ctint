// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
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
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t &qmc_config;

    // Container for the accumulation
    gf_view<mesh::dlr_imtime, tensor_valued<1>> chiAB_tau_;

    // Case types for operator block structure
    enum class chi_case_t { AAAA, AABB, ABAB };

    // A single entry in a chi group
    struct chi_entry_t {
      long pair_idx;           // target index in chiAB_tau_
      dcomplex coef;           // coef_A * coef_B (ABAB minus sign absorbed)
      int idx_cdag_A, idx_c_A; // A-side orbital indices
      int idx_cdag_B, idx_c_B; // B-side orbital indices
    };

    // A group of entries sharing the same case type and block indices
    struct chi_group_t {
      chi_case_t case_type;
      int bl_det1, bl_det2;              // block indices for determinant(s)
      std::vector<chi_entry_t> entries;
      // Pre-allocated scratch arrays of shape (L, entries.size()), filled each MC step
      nda::array<c_t, 2> c_A, c_B;
      nda::array<cdag_t, 2> cdag_A, cdag_B;
    };

    // Grouped operator pairs
    std::vector<chi_group_t> groups_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Precomputed tau points from the DLR imtime mesh
    long L_;
    std::vector<tau_t> tau_points_;
  };

} // namespace triqs_ctint::measures
