// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./chiAB_tau.hpp"
#include "./../types.hpp"

using namespace triqs::utility;
using triqs::operators::bilinear_type;

namespace triqs_ctint::measures {

  chiAB_tau::chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_) {

    if (params.chi_ops.empty()) TRIQS_RUNTIME_ERROR << " Empty operator pair list detected in chiAB measurement \n";

    // Parse operator pairs
    using op_term_t = std::tuple<dcomplex, bilinear_type, std::pair<int, int>, std::pair<int, int>>;
    std::vector<std::pair<std::vector<op_term_t>, std::vector<op_term_t>>> op_pairs;
    for (auto const &[A, B] : params.chi_ops) op_pairs.emplace_back(get_terms(A, params.gf_struct), get_terms(B, params.gf_struct));

    // Group all (target_idx, A_term, B_term) triples by (op_type, case_type, block_indices)
    using group_key_t = std::tuple<quartic_op_type, quartic_case_t, int, int>;
    std::map<group_key_t, quartic_group_t> group_map;

    for (auto [target_idx, pair] : itertools::enumerate(op_pairs)) {
      auto &[A, B] = pair;
      for (auto &[coef_B, type_B, bl_pair_B, idx_pair_B] : B) {
        auto [bl1_B, bl2_B] = bl_pair_B;
        auto [idx1_B, idx2_B] = idx_pair_B;

        for (auto &[coef_A, type_A, bl_pair_A, idx_pair_A] : A) {
          auto [bl1_A, bl2_A] = bl_pair_A;
          auto [idx1_A, idx2_A] = idx_pair_A;

          // Determine the quartic operator type from the bilinear types.
          // Valid: (c†c, c†c), (c†c†, cc), (cc, c†c†). Others have odd fermion number and are skipped.
          quartic_op_type op_type;
          // Block indices for classify: (row1_bl, col1_bl, row2_bl, col2_bl)
          int bl_row1, bl_col1, bl_row2, bl_col2;
          // Orbital indices for quartic_entry_t: (cdag_A=row1, c_A=col1, cdag_B=row2, c_B=col2)
          int e_idx_cdag_A, e_idx_c_A, e_idx_cdag_B, e_idx_c_B;

          if (type_A == bilinear_type::cdag_c && type_B == bilinear_type::cdag_c) {
            // Normal: A=(c†,c), B=(c†,c) — rows from c†'s, cols from c's
            op_type    = quartic_op_type::normal;
            bl_row1    = bl1_A;  bl_col1 = bl2_A;   // A's c† and c
            bl_row2    = bl1_B;  bl_col2 = bl2_B;   // B's c† and c
            e_idx_cdag_A = idx1_A; e_idx_c_A = idx2_A;
            e_idx_cdag_B = idx1_B; e_idx_c_B = idx2_B;
          } else if (type_A == bilinear_type::cdag_cdag && type_B == bilinear_type::c_c) {
            // Anomalous: A=(c†,c†), B=(c,c) — both rows from A, both cols from B
            op_type    = quartic_op_type::anom_cdcd_cc;
            bl_row1    = bl1_A;  bl_col1 = bl1_B;   // A's first c† paired with B's first c
            bl_row2    = bl2_A;  bl_col2 = bl2_B;   // A's second c† paired with B's second c
            e_idx_cdag_A = idx1_A; e_idx_c_A = idx1_B;
            e_idx_cdag_B = idx2_A; e_idx_c_B = idx2_B;
          } else if (type_A == bilinear_type::c_c && type_B == bilinear_type::cdag_cdag) {
            // Anomalous reversed: A=(c,c), B=(c†,c†) — both cols from A, both rows from B
            op_type    = quartic_op_type::anom_cc_cdcd;
            bl_row1    = bl1_B;  bl_col1 = bl1_A;   // B's first c† paired with A's first c
            bl_row2    = bl2_B;  bl_col2 = bl2_A;   // B's second c† paired with A's second c
            e_idx_cdag_A = idx1_B; e_idx_c_A = idx1_A;
            e_idx_cdag_B = idx2_B; e_idx_c_B = idx2_A;
          } else {
            continue; // Invalid combination (odd fermion number) — skip
          }

          // Anomalous ordering (c†c†cc) picks up an extra (-1) relative to normal ordering (c†cc†c)
          // because the Wick contraction of c†c†cc relates to the determinant with an extra sign:
          // det(M) equals the Wick contraction for c†cc†c, but -det(M) for c†c†cc.
          auto combined_coef = coef_A * coef_B * (op_type == quartic_op_type::normal ? 1.0 : -1.0);
          auto info          = classify_quartic_blocks(bl_row1, bl_col1, bl_row2, bl_col2, combined_coef);
          if (!info) continue;

          group_key_t key{op_type, info->case_type, info->bl_det1, info->bl_det2};
          auto &grp     = group_map[key];
          grp.case_type = info->case_type;
          grp.op_type   = op_type;
          grp.bl_det1   = info->bl_det1;
          grp.bl_det2   = info->bl_det2;
          grp.entries.push_back({static_cast<long>(target_idx), info->coef, e_idx_cdag_A, e_idx_c_A, e_idx_cdag_B, e_idx_c_B});
        }
      }
    }

    // Init measurement container and capture view
    mesh::dlr_imtime tau_mesh{params_.beta, Boson, params_.dlr_wmax, params_.dlr_eps, true};
    results->chiAB_tau = gf<mesh::dlr_imtime, tensor_valued<1>>{tau_mesh, make_shape(op_pairs.size())};
    chiAB_tau_.rebind(results->chiAB_tau.value());
    chiAB_tau_() = 0;

    // Precompute tau points from the DLR mesh
    L_ = tau_mesh.size();
    tau_points_.reserve(L_);
    for (auto tau : tau_mesh) tau_points_.push_back(make_tau_t(double(tau)));

    // Flatten map to vector, pre-allocate scratch arrays, and fill constant side
    groups_.reserve(group_map.size());
    for (auto &[key, grp] : group_map) {
      long E = static_cast<long>(grp.entries.size());
      grp.c_A.resize(L_, E);
      grp.cdag_A.resize(L_, E);
      grp.c_B.resize(L_, E);
      grp.cdag_B.resize(L_, E);

      // Fill the tau-independent (constant) side of the operator arrays
      for (long l = 0; l < L_; ++l) {
        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          switch (grp.op_type) {
            case quartic_op_type::normal:
              // B-side (time 0): one row + one col
              grp.c_B(l, e)    = c_t{tau_t::get_zero(), entry.idx_c_B};
              grp.cdag_B(l, e) = cdag_t{tau_t::get_zero_plus(), entry.idx_cdag_B};
              break;
            case quartic_op_type::anom_cdcd_cc:
              // B-side (time 0): two cols
              grp.c_A(l, e) = c_t{tau_t::get_zero_plus(), entry.idx_c_A};
              grp.c_B(l, e) = c_t{tau_t::get_zero(), entry.idx_c_B};
              break;
            case quartic_op_type::anom_cc_cdcd:
              // B-side (time 0): two rows
              grp.cdag_A(l, e) = cdag_t{tau_t::get_zero_plus(), entry.idx_cdag_A};
              grp.cdag_B(l, e) = cdag_t{tau_t::get_zero(), entry.idx_cdag_B};
              break;
          }
        }
      }
      groups_.push_back(std::move(grp));
    }
  }

  void chiAB_tau::accumulate(mc_weight_t sign) {
    Z += sign;

    for (auto &grp : groups_) {
      long E = static_cast<long>(grp.entries.size());

      // Fill tau-varying side of the operator arrays
      for (long l = 0; l < L_; ++l) {
        auto tau     = tau_points_[l];
        auto tau_lo  = tau_t{tau.n + 2}; // second operator in A (earlier infinitesimal)
        auto tau_hi  = tau_t{tau.n + 3}; // first operator in A (later infinitesimal)
        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          switch (grp.op_type) {
            case quartic_op_type::normal:
              // A-side (time tau): one row + one col
              grp.cdag_A(l, e) = cdag_t{tau_hi, entry.idx_cdag_A};
              grp.c_A(l, e)    = c_t{tau_lo, entry.idx_c_A};
              break;
            case quartic_op_type::anom_cdcd_cc:
              // A-side (time tau): two rows
              grp.cdag_A(l, e) = cdag_t{tau_hi, entry.idx_cdag_A};
              grp.cdag_B(l, e) = cdag_t{tau_lo, entry.idx_cdag_B};
              break;
            case quartic_op_type::anom_cc_cdcd:
              // A-side (time tau): two cols
              grp.c_A(l, e) = c_t{tau_hi, entry.idx_c_A};
              grp.c_B(l, e) = c_t{tau_lo, entry.idx_c_B};
              break;
          }
        }
      }

      // Scatter ratios into chiAB_tau_ (l-outer for row-major cache locality)
      auto scatter = [&](auto const &ratios) {
        for (long l = 0; l < L_; ++l)
          for (long e = 0; e < E; ++e) chiAB_tau_.data()(l, grp.entries[e].target_idx) += grp.entries[e].coef * sign * ratios(l, e);
      };

      switch (grp.case_type) {
        case quartic_case_t::AAAA: {
          scatter(qmc_config.dets[grp.bl_det1].insert2_ratios(0, 1, 0, 1, grp.c_A, grp.c_B, grp.cdag_A, grp.cdag_B));
          break;
        }
        case quartic_case_t::AABB: {
          auto r1 = qmc_config.dets[grp.bl_det1].insert_ratios(0, 0, grp.c_A, grp.cdag_A);
          auto r2 = qmc_config.dets[grp.bl_det2].insert_ratios(0, 0, grp.c_B, grp.cdag_B);
          r1 *= r2;
          scatter(r1);
          break;
        }
        case quartic_case_t::ABAB: {
          // det_1: (c_B, cdag_A), det_2: (c_A, cdag_B) — cross-pair into each det
          auto r1 = qmc_config.dets[grp.bl_det1].insert_ratios(0, 0, grp.c_B, grp.cdag_A);
          auto r2 = qmc_config.dets[grp.bl_det2].insert_ratios(0, 0, grp.c_A, grp.cdag_B);
          r1 *= r2;
          scatter(r1);
          break;
        }
      }
    }
  }

  void chiAB_tau::collect_results(mpi::communicator const &comm) {
    Z          = mpi::all_reduce(Z, comm);
    chiAB_tau_ = mpi::all_reduce(chiAB_tau_, comm);
    chiAB_tau_ = chiAB_tau_ / Z;
  }

} // namespace triqs_ctint::measures
