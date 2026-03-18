// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/stat/lin_binning.hpp>
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of static expectation values $\langle C_i \rangle$ by operator insertion with tau-averaging.
  * Each C_i is a many_body_operator containing constant, bilinear (c†c), and quartic (c†c†cc) terms.
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

    // Case types for quartic operator block structure
    enum class case_t { AAAA, AABB, ABAB };

    // Bilinear entry: single insertion into one det
    struct bilinear_entry_t {
      long obs_idx;
      dcomplex coef;
      int idx_cdag, idx_c;
    };
    struct bilinear_group_t {
      int bl_det;
      std::vector<bilinear_entry_t> entries;
    };

    // Quartic entry: double insertion or product of insertions
    struct quartic_entry_t {
      long obs_idx;
      dcomplex coef;
      int idx_cdag_A, idx_c_A; // pair A: (m[0], m[3])
      int idx_cdag_B, idx_c_B; // pair B: (m[1], m[2])
    };
    struct quartic_group_t {
      case_t case_type;
      int bl_det1, bl_det2;
      std::vector<quartic_entry_t> entries;
    };

    // Grouped operator terms
    std::vector<bilinear_group_t> bilinear_groups_;
    std::vector<quartic_group_t> quartic_groups_;

    // Constant contributions (degree-0 monomials)
    nda::array<dcomplex, 1> constant_parts_;

    // Tau averaging
    long L_;                         // number of tau points
    std::vector<tau_t> tau_points_;  // precomputed uniform tau grid

    mc_weight_t Z = 0.0;
    long N_        = 0; // step counter

    // Per-step contributions and linear binning for error analysis
    nda::array<dcomplex, 1> step_contrib_;
    std::vector<triqs::stat::lin_binning<dcomplex>> obs_bins_;
  };

} // namespace triqs_ctint::measures
