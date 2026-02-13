// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/mesh/matsubara_freq.hpp>
#include <tuple>
#include <vector>

namespace triqs_ctint::measures {

  /// Build target Matsubara frequency vector and index map from a set of unique frequencies.
  /// Returns:
  ///   target_mf: vector of Matsubara frequencies for NFFT type 3
  ///   idx_offset: negated Matsubara index of first unique frequency
  ///   to_idx: flat lookup mapping (n + idx_offset) -> position in target_mf
  template <typename Set>
  inline auto build_index_map(Set const &unique_w_set) {
    std::vector<triqs::mesh::matsubara_freq> target_mf(unique_w_set.begin(), unique_w_set.end());

    long idx_offset = -target_mf.front().n;
    std::vector<long> to_idx(target_mf.back().n - target_mf.front().n + 1, -1);

    for (size_t k = 0; k < target_mf.size(); ++k) to_idx[target_mf[k].n + idx_offset] = static_cast<long>(k);

    return std::tuple{std::move(target_mf), idx_offset, std::move(to_idx)};
  }

} // namespace triqs_ctint::measures
