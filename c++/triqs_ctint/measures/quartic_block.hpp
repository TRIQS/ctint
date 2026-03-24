// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include <optional>

namespace triqs_ctint::measures {

  // Quartic operator block-structure classification
  enum class quartic_case_t { AAAA, AABB, ABAB };

  // Classification of bilinear operator pair type for quartic measurements
  enum class quartic_op_type { normal, anom_cdcd_cc, anom_cc_cdcd };

  struct quartic_block_info_t {
    quartic_case_t case_type;
    int bl_det1, bl_det2;
    dcomplex coef; // ABAB sign absorbed
  };

  // Classify a quartic c†c c†c term by block indices. Returns nullopt if the combination vanishes.
  inline std::optional<quartic_block_info_t>
  classify_quartic_blocks(int bl_cdag_A, int bl_c_A, int bl_cdag_B, int bl_c_B, dcomplex coef) {
    bool is_AABB = (bl_cdag_A == bl_c_A && bl_cdag_B == bl_c_B);
    bool is_ABAB = (bl_cdag_A == bl_c_B && bl_cdag_B == bl_c_A);
    if (is_AABB && is_ABAB) return quartic_block_info_t{quartic_case_t::AAAA, bl_cdag_A, bl_cdag_A, coef};
    if (is_AABB)            return quartic_block_info_t{quartic_case_t::AABB, bl_cdag_A, bl_cdag_B, coef};
    if (is_ABAB)            return quartic_block_info_t{quartic_case_t::ABAB, bl_cdag_A, bl_c_A, -coef};
    return std::nullopt;
  }

  // Quartic entry: orbital indices for one c†c c†c term, targeting a specific output slot
  struct quartic_entry_t {
    long target_idx; // output index (pair_idx for chiAB, obs_idx for static_obs)
    dcomplex coef;
    int idx_cdag_A, idx_c_A;
    int idx_cdag_B, idx_c_B;
  };

  // Group of quartic entries sharing case type, block indices, and operator type
  struct quartic_group_t {
    quartic_case_t case_type;
    quartic_op_type op_type = quartic_op_type::normal;
    int bl_det1, bl_det2;
    std::vector<quartic_entry_t> entries;
    // Pre-allocated scratch arrays of shape (L, entries.size()), filled each MC step
    nda::array<c_t, 2> c_A, c_B;
    nda::array<cdag_t, 2> cdag_A, cdag_B;
  };

} // namespace triqs_ctint::measures
