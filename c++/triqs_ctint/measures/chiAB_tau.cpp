// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./chiAB_tau.hpp"
#include "./../types.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  chiAB_tau::chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), tau_mesh{params_.beta, Boson, params_.dlr_wmax, params_.dlr_eps} {

    if (params.chi_ops.empty()) TRIQS_RUNTIME_ERROR << " Empty operator pair list detected in chiAB measurement \n";

    // Parse operator pairs
    using op_term_t = std::tuple<dcomplex, std::pair<int, int>, std::pair<int, int>>;
    std::vector<std::pair<std::vector<op_term_t>, std::vector<op_term_t>>> op_pairs;
    for (auto const &[A, B] : params.chi_ops) op_pairs.emplace_back(get_terms(A, params.gf_struct), get_terms(B, params.gf_struct));

    // Group all (pair_idx, A_term, B_term) triples by (case_type, block_indices)
    // Use a map keyed by (case_type, bl_det1, bl_det2)
    using group_key_t = std::tuple<chi_case_t, int, int>;
    std::map<group_key_t, chi_group_t> group_map;

    for (auto [pair_idx, pair] : itertools::enumerate(op_pairs)) {
      auto &[A, B] = pair;
      for (auto &[coef_B, bl_pair_B, idx_pair_B] : B) {
        auto [bl_cdag_B, bl_c_B] = bl_pair_B;
        auto [idx_cdag_B, idx_c_B] = idx_pair_B;

        for (auto &[coef_A, bl_pair_A, idx_pair_A] : A) {
          auto [bl_cdag_A, bl_c_A] = bl_pair_A;
          auto [idx_cdag_A, idx_c_A] = idx_pair_A;

          bool is_AABB = (bl_cdag_A == bl_c_A && bl_cdag_B == bl_c_B);
          bool is_ABAB = (bl_cdag_A == bl_c_B && bl_cdag_B == bl_c_A);
          bool is_AAAA = is_AABB && is_ABAB;

          chi_case_t case_type;
          int bl_det1, bl_det2;
          dcomplex coef = coef_A * coef_B;

          if (is_AAAA) {
            case_type = chi_case_t::AAAA;
            bl_det1 = bl_cdag_A;
            bl_det2 = bl_cdag_A; // unused, same det
          } else if (is_AABB) {
            case_type = chi_case_t::AABB;
            bl_det1 = bl_cdag_A;
            bl_det2 = bl_cdag_B;
          } else if (is_ABAB) {
            case_type = chi_case_t::ABAB;
            bl_det1 = bl_cdag_A;
            bl_det2 = bl_c_A;
            coef = -coef; // absorb minus sign from operator swap
          } else {
            continue; // other block combinations vanish
          }

          group_key_t key{case_type, bl_det1, bl_det2};
          auto &grp = group_map[key];
          grp.case_type = case_type;
          grp.bl_det1 = bl_det1;
          grp.bl_det2 = bl_det2;
          grp.entries.push_back({static_cast<long>(pair_idx), coef, idx_cdag_A, idx_c_A, idx_cdag_B, idx_c_B});
        }
      }
    }

    // Flatten map to vector
    groups_.reserve(group_map.size());
    for (auto &[key, grp] : group_map) groups_.push_back(std::move(grp));

    // Init measurement container and capture view
    results->chiAB_tau = gf<mesh::dlr_imtime, tensor_valued<1>>{tau_mesh, make_shape(op_pairs.size())};
    chiAB_tau_.rebind(results->chiAB_tau.value());
    chiAB_tau_() = 0;
  }

  void chiAB_tau::accumulate(mc_weight_t sign) {
    Z += sign;

    long L = tau_mesh.size();

    for (auto const &grp : groups_) {
      long E = static_cast<long>(grp.entries.size());
      auto &det_1 = qmc_config.dets[grp.bl_det1];

      if (grp.case_type == chi_case_t::AAAA) {
        // Build rank-2 arrays (L, E) for A-side (tau varies), rank-1 (E) for B-side (tau=0 fixed)
        nda::array<c_t, 2> c_A(L, E);
        nda::array<cdag_t, 2> cdag_A(L, E);
        nda::array<c_t, 1> c_B(E);
        nda::array<cdag_t, 1> cdag_B(E);

        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          c_B(e) = c_t{tau_t::get_zero(), entry.idx_c_B};
          cdag_B(e) = cdag_t{tau_t::get_zero_plus(), entry.idx_cdag_B};
        }

        long l = 0;
        for (auto tau : tau_mesh) {
          auto tau_point = make_tau_t(double(tau));
          auto taup_point = tau_point;
          tau_point.n += 3;
          taup_point.n += 2;
          for (long e = 0; e < E; ++e) {
            auto const &entry = grp.entries[e];
            c_A(l, e) = c_t{taup_point, entry.idx_c_A};
            cdag_A(l, e) = cdag_t{tau_point, entry.idx_cdag_A};
          }
          ++l;
        }

        // Single batched call with broadcast: (L,E) × (E) → (L,E)
        auto ratios = det_1.insert2_ratios(0, 1, 0, 1, c_A, c_B, cdag_A, cdag_B);

        // Scatter results
        l = 0;
        for (auto tau : tau_mesh) {
          for (long e = 0; e < E; ++e) { chiAB_tau_[tau](grp.entries[e].pair_idx) += grp.entries[e].coef * sign * ratios(l, e); }
          ++l;
        }

      } else if (grp.case_type == chi_case_t::AABB) {
        auto &det_2 = qmc_config.dets[grp.bl_det2];

        // A-side: rank-2 (L, E), B-side: rank-1 (E)
        nda::array<c_t, 2> c_A(L, E);
        nda::array<cdag_t, 2> cdag_A(L, E);
        nda::array<c_t, 1> c_B(E);
        nda::array<cdag_t, 1> cdag_B(E);

        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          c_B(e) = c_t{tau_t::get_zero(), entry.idx_c_B};
          cdag_B(e) = cdag_t{tau_t::get_zero_plus(), entry.idx_cdag_B};
        }

        long l = 0;
        for (auto tau : tau_mesh) {
          auto tau_point = make_tau_t(double(tau));
          auto taup_point = tau_point;
          tau_point.n += 3;
          taup_point.n += 2;
          for (long e = 0; e < E; ++e) {
            auto const &entry = grp.entries[e];
            c_A(l, e) = c_t{taup_point, entry.idx_c_A};
            cdag_A(l, e) = cdag_t{tau_point, entry.idx_cdag_A};
          }
          ++l;
        }

        // det_1: A-side, rank-2 (L, E) × rank-2 (L, E)
        auto ratios_1 = det_1.insert_ratios(0, 0, c_A, cdag_A);
        // det_2: B-side, rank-1 (E) × rank-1 (E)
        auto ratios_2 = det_2.insert_ratios(0, 0, c_B, cdag_B);

        // Scatter with broadcast: ratios_1(l,e) * ratios_2(e)
        l = 0;
        for (auto tau : tau_mesh) {
          for (long e = 0; e < E; ++e) {
            chiAB_tau_[tau](grp.entries[e].pair_idx) += grp.entries[e].coef * sign * ratios_1(l, e) * ratios_2(e);
          }
          ++l;
        }

      } else { // ABAB
        auto &det_2 = qmc_config.dets[grp.bl_det2];

        // ABAB: det_1 gets (c_B, cdag_A), det_2 gets (c_A, cdag_B)
        // c_B has fixed tau=0, cdag_A varies with tau
        // c_A varies with tau, cdag_B has fixed tau=0+
        // Both det calls need rank-2 arrays (replicate fixed-tau side)
        nda::array<c_t, 2> c_B_rep(L, E);    // replicated across L
        nda::array<cdag_t, 2> cdag_A(L, E);   // varies with tau
        nda::array<c_t, 2> c_A(L, E);         // varies with tau
        nda::array<cdag_t, 2> cdag_B_rep(L, E); // replicated across L

        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          auto c_B_e = c_t{tau_t::get_zero(), entry.idx_c_B};
          auto cdag_B_e = cdag_t{tau_t::get_zero_plus(), entry.idx_cdag_B};
          for (long l2 = 0; l2 < L; ++l2) {
            c_B_rep(l2, e) = c_B_e;
            cdag_B_rep(l2, e) = cdag_B_e;
          }
        }

        long l = 0;
        for (auto tau : tau_mesh) {
          auto tau_point = make_tau_t(double(tau));
          auto taup_point = tau_point;
          tau_point.n += 3;
          taup_point.n += 2;
          for (long e = 0; e < E; ++e) {
            auto const &entry = grp.entries[e];
            c_A(l, e) = c_t{taup_point, entry.idx_c_A};
            cdag_A(l, e) = cdag_t{tau_point, entry.idx_cdag_A};
          }
          ++l;
        }

        // det_1: (c_B, cdag_A) — swapped operators
        auto ratios_1 = det_1.insert_ratios(0, 0, c_B_rep, cdag_A);
        // det_2: (c_A, cdag_B)
        auto ratios_2 = det_2.insert_ratios(0, 0, c_A, cdag_B_rep);

        // Scatter (minus sign already absorbed into coef)
        l = 0;
        for (auto tau : tau_mesh) {
          for (long e = 0; e < E; ++e) {
            chiAB_tau_[tau](grp.entries[e].pair_idx) += grp.entries[e].coef * sign * ratios_1(l, e) * ratios_2(l, e);
          }
          ++l;
        }
      }
    }
  }

  void chiAB_tau::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z          = mpi::all_reduce(Z, comm);
    chiAB_tau_ = mpi::all_reduce(chiAB_tau_, comm);
    chiAB_tau_ = chiAB_tau_ / Z;
  }

} // namespace triqs_ctint::measures
