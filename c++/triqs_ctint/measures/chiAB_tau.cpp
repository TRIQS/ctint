// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./chiAB_tau.hpp"
#include "./../types.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  chiAB_tau::chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_) {

    if (params.chi_ops.empty()) TRIQS_RUNTIME_ERROR << " Empty operator pair list detected in chiAB measurement \n";

    // Parse operator pairs
    using op_term_t = std::tuple<dcomplex, std::pair<int, int>, std::pair<int, int>>;
    std::vector<std::pair<std::vector<op_term_t>, std::vector<op_term_t>>> op_pairs;
    for (auto const &[A, B] : params.chi_ops) op_pairs.emplace_back(get_terms(A, params.gf_struct), get_terms(B, params.gf_struct));

    // Group all (target_idx, A_term, B_term) triples by (case_type, block_indices)
    using group_key_t = std::tuple<quartic_case_t, int, int>;
    std::map<group_key_t, quartic_group_t> group_map;

    for (auto [target_idx, pair] : itertools::enumerate(op_pairs)) {
      auto &[A, B] = pair;
      for (auto &[coef_B, bl_pair_B, idx_pair_B] : B) {
        auto [bl_cdag_B, bl_c_B] = bl_pair_B;
        auto [idx_cdag_B, idx_c_B] = idx_pair_B;

        for (auto &[coef_A, bl_pair_A, idx_pair_A] : A) {
          auto [bl_cdag_A, bl_c_A] = bl_pair_A;
          auto [idx_cdag_A, idx_c_A] = idx_pair_A;

          auto info = classify_quartic_blocks(bl_cdag_A, bl_c_A, bl_cdag_B, bl_c_B, coef_A * coef_B);
          if (!info) continue;

          group_key_t key{info->case_type, info->bl_det1, info->bl_det2};
          auto &grp     = group_map[key];
          grp.case_type = info->case_type;
          grp.bl_det1   = info->bl_det1;
          grp.bl_det2   = info->bl_det2;
          grp.entries.push_back({static_cast<long>(target_idx), info->coef, idx_cdag_A, idx_c_A, idx_cdag_B, idx_c_B});
        }
      }
    }

    // Init measurement container and capture view
    mesh::dlr_imtime tau_mesh{params_.beta, Boson, params_.dlr_wmax, params_.dlr_eps};
    results->chiAB_tau = gf<mesh::dlr_imtime, tensor_valued<1>>{tau_mesh, make_shape(op_pairs.size())};
    chiAB_tau_.rebind(results->chiAB_tau.value());
    chiAB_tau_() = 0;

    // Precompute tau points from the DLR mesh
    L_ = tau_mesh.size();
    tau_points_.reserve(L_);
    for (auto tau : tau_mesh) tau_points_.push_back(make_tau_t(double(tau)));

    // Flatten map to vector, pre-allocate scratch arrays, and fill constant B-side
    groups_.reserve(group_map.size());
    for (auto &[key, grp] : group_map) {
      long E = static_cast<long>(grp.entries.size());
      grp.c_A.resize(L_, E);
      grp.cdag_A.resize(L_, E);
      grp.c_B.resize(L_, E);
      grp.cdag_B.resize(L_, E);
      // B-side operators are at tau=0 and never change — fill once
      for (long l = 0; l < L_; ++l) {
        for (long e = 0; e < E; ++e) {
          grp.c_B(l, e)    = c_t{tau_t::get_zero(), grp.entries[e].idx_c_B};
          grp.cdag_B(l, e) = cdag_t{tau_t::get_zero_plus(), grp.entries[e].idx_cdag_B};
        }
      }
      groups_.push_back(std::move(grp));
    }
  }

  void chiAB_tau::accumulate(mc_weight_t sign) {
    Z += sign;

    for (auto &grp : groups_) {
      long E = static_cast<long>(grp.entries.size());

      // Fill A-side scratch arrays (vary with tau; B-side is constant, filled in constructor)
      for (long l = 0; l < L_; ++l) {
        auto tau      = tau_points_[l];
        auto tau_cA   = tau_t{tau.n + 2};
        auto tau_cdA  = tau_t{tau.n + 3};
        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          grp.c_A(l, e)    = c_t{tau_cA, entry.idx_c_A};
          grp.cdag_A(l, e) = cdag_t{tau_cdA, entry.idx_cdag_A};
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
