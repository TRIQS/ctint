// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <nda/nda.hpp>
#include <triqs/mesh/matsubara_freq.hpp>
#include <vector>

namespace triqs_ctint::measures {

  /// Build target Matsubara frequency array and index map from sorted unique frequencies.
  /// target_mf: shape (1, n_unique) — Matsubara frequencies for NFFT type 3
  /// idx_offset: negated Matsubara index of first unique frequency
  /// to_idx: flat lookup mapping (n + idx_offset) -> position in unique_w
  inline void build_index_map(std::vector<triqs::mesh::matsubara_freq> const &unique_w,
                              nda::array<triqs::mesh::matsubara_freq, 2> &target_mf, long &idx_offset, std::vector<long> &to_idx) {
    auto n_un = static_cast<int64_t>(unique_w.size());
    target_mf.resize(1, n_un);
    idx_offset = -unique_w.front().n;
    to_idx.assign(unique_w.back().n - unique_w.front().n + 1, -1);
    for (int64_t k = 0; k < n_un; ++k) {
      target_mf(0, k)                    = unique_w[k];
      to_idx[unique_w[k].n + idx_offset] = k;
    }
  }

} // namespace triqs_ctint::measures
