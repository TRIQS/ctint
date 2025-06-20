// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./chiAB_tau.hpp"
#include "./../types.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  chiAB_tau::chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), tau_mesh{params_.beta, Boson, params_.n_tau_chi2} {

    if (params.chi_A_vec.size() == 0 or params.chi_B_vec.size() == 0)
      TRIQS_RUNTIME_ERROR << " Empty operator vector detected in chiAB measurement \n";

    for (auto A : params.chi_A_vec) A_vec.emplace_back(get_terms(A, params.gf_struct));
    for (auto B : params.chi_B_vec) B_vec.emplace_back(get_terms(B, params.gf_struct));

    // Init measurement container and capture view
    results->chiAB_tau = gf<imtime>{tau_mesh, make_shape(A_vec.size(), B_vec.size())};
    chiAB_tau_.rebind(results->chiAB_tau.value());
    chiAB_tau_() = 0;
  }

  void chiAB_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    for (auto [j, B] : enumerate(B_vec))
      for (auto &[coef_B, bl_pair_B, idx_pair_B] : B) {

        auto [idx_cdag_B, idx_c_B] = idx_pair_B;
        auto cdag_B                = cdag_t{tau_t::get_zero_plus(), idx_cdag_B};
        auto c_B                   = c_t{tau_t::get_zero(), idx_c_B};
        auto [bl_cdag_B, bl_c_B]   = bl_pair_B;

        for (auto [i, A] : enumerate(A_vec))
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

            for (auto tau : tau_mesh) {
              auto tau_point  = make_tau_t(double(tau));
              auto taup_point = tau_point;

              // Time-Ordering
              tau_point.n += 3;  // tau -> tau^{+++}
              taup_point.n += 2; // taup -> tau^{++}

              if (tau.index() == params.n_tau_chi2 - 1) {
                tau_point  = tau_t::get_beta_minus();
                taup_point = tau_t::get_beta_minus_minus();
              }

              cdag_A.tau = tau_point;
              c_A.tau    = taup_point;

              // TODO Extend for measurement of pairing response functions
              if (is_AAAA)
                chiAB_tau_[tau](i, j) += coef_A * coef_B * sign * det_1.try_insert2(0, 1, 0, 1, c_A, c_B, cdag_A, cdag_B);
              else if (is_AABB)
                chiAB_tau_[tau](i, j) += coef_A * coef_B * sign * det_1.try_insert(0, 0, c_A, cdag_A) * det_2.try_insert(0, 0, c_B, cdag_B);
              else if (is_ABAB) // In this case we have to swap c_A and c_B, which leads to a minus sign
                chiAB_tau_[tau](i, j) -= coef_A * coef_B * sign * det_1.try_insert(0, 0, c_B, cdag_A) * det_2.try_insert(0, 0, c_A, cdag_B);
              // Note: All other block combinations vanish

              det_1.reject_last_try();
              det_2.reject_last_try();
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
