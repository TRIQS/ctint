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

    for (auto const &[A, B] : params.chi_ops) op_pairs.emplace_back(get_terms(A, params.gf_struct), get_terms(B, params.gf_struct));

    // Init measurement container and capture view
    results->chiAB_tau = gf<mesh::dlr_imtime, tensor_valued<1>>{tau_mesh, make_shape(op_pairs.size())};
    chiAB_tau_.rebind(results->chiAB_tau.value());
    chiAB_tau_() = 0;
  }

  void chiAB_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    long K = tau_mesh.size();

    for (auto [i, pair] : enumerate(op_pairs)) {
      auto &[A, B] = pair;
      for (auto &[coef_B, bl_pair_B, idx_pair_B] : B) {

        auto [idx_cdag_B, idx_c_B] = idx_pair_B;
        auto cdag_B                = cdag_t{tau_t::get_zero_plus(), idx_cdag_B};
        auto c_B                   = c_t{tau_t::get_zero(), idx_c_B};
        auto [bl_cdag_B, bl_c_B]   = bl_pair_B;

        for (auto &[coef_A, bl_pair_A, idx_pair_A] : A) {

          auto [idx_cdag_A, idx_c_A] = idx_pair_A;
          auto cdag_A                = cdag_t{tau_t::get_zero_plus(), idx_cdag_A};
          auto c_A                   = c_t{tau_t::get_zero(), idx_c_A};
          auto [bl_cdag_A, bl_c_A]   = bl_pair_A;

          // Determine the block order in the quartet of operators
          bool is_AABB = (bl_cdag_A == bl_c_A && bl_cdag_B == bl_c_B);
          bool is_ABAB = (bl_cdag_A == bl_c_B && bl_cdag_B == bl_c_A);
          bool is_AAAA = is_AABB && is_ABAB;

          // Define the determinants to operate on
          auto &det_1 = qmc_config.dets[bl_cdag_A];
          auto &det_2 = is_AABB ? qmc_config.dets[bl_cdag_B] : qmc_config.dets[bl_c_A];

          // Build vectors of operator arguments for all tau points (batch insertion)
          std::vector<c_t> c_A_vec(K), c_B_vec(K);
          std::vector<cdag_t> cdag_A_vec(K), cdag_B_vec(K);

          long k = 0;
          for (auto tau : tau_mesh) {
            auto tau_point  = make_tau_t(double(tau));
            auto taup_point = tau_point;
            tau_point.n += 3;  // tau -> tau^{+++}
            taup_point.n += 2; // taup -> tau^{++}

            c_A_vec[k]    = c_A;
            c_A_vec[k].tau = taup_point;
            cdag_A_vec[k] = cdag_A;
            cdag_A_vec[k].tau = tau_point;
            c_B_vec[k]    = c_B;
            cdag_B_vec[k] = cdag_B;
            ++k;
          }

          // TODO Extend for measurement of pairing response functions
          if (is_AAAA) {
            auto ratios = det_1.insert2_ratios(0, 1, 0, 1, c_A_vec, c_B_vec, cdag_A_vec, cdag_B_vec);
            k = 0;
            for (auto tau : tau_mesh) { chiAB_tau_[tau](i) += coef_A * coef_B * sign * ratios(k++); }
          } else if (is_AABB) {
            auto ratios_1 = det_1.insert_ratios(0, 0, c_A_vec, cdag_A_vec);
            auto ratios_2 = det_2.insert_ratios(0, 0, c_B_vec, cdag_B_vec);
            k = 0;
            for (auto tau : tau_mesh) { chiAB_tau_[tau](i) += coef_A * coef_B * sign * ratios_1(k) * ratios_2(k); ++k; }
          } else if (is_ABAB) {
            // Swap c_A and c_B which leads to a minus sign
            auto ratios_1 = det_1.insert_ratios(0, 0, c_B_vec, cdag_A_vec);
            auto ratios_2 = det_2.insert_ratios(0, 0, c_A_vec, cdag_B_vec);
            k = 0;
            for (auto tau : tau_mesh) { chiAB_tau_[tau](i) -= coef_A * coef_B * sign * ratios_1(k) * ratios_2(k); ++k; }
          }
          // Note: All other block combinations vanish
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
