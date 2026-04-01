// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./chiAB_tau.hpp"

namespace triqs_ctint::measures {

  chiAB_tau::chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), n_blocks_(static_cast<int>(params_.gf_struct.size())) {

    if (params.chi_ops.empty()) TRIQS_RUNTIME_ERROR << " Empty operator pair list detected in chiAB measurement \n";

    // Group all (target_idx, A_term, B_term) triples by block signature
    std::map<block_signature_t, term_group_t> group_map;

    for (auto [target_idx, pair] : itertools::enumerate(params.chi_ops)) {
      auto &[A_op, B_op] = pair;

      for (auto const &term_a : A_op) {
        for (auto const &term_b : B_op) {
          auto const &ma = term_a.monomial;
          auto const &mb = term_b.monomial;

          // Check total fermion balance
          int n_cdag = 0, n_c = 0;
          for (auto const &op : ma) op.dagger ? ++n_cdag : ++n_c;
          for (auto const &op : mb) op.dagger ? ++n_cdag : ++n_c;
          if (n_cdag != n_c) continue;

          // Skip constant terms
          if (ma.empty() || mb.empty()) continue;

          auto result =
             make_chiAB_term(ma, mb, dcomplex(term_a.coef) * dcomplex(term_b.coef), static_cast<long>(target_idx), params.gf_struct, n_blocks_);
          if (!result) continue;
          auto &[oterm, sig]       = *result;
          group_map[sig].signature = sig;
          group_map[sig].terms.push_back(std::move(oterm));
        }
      }
    }

    // Init measurement container and capture view
    mesh::dlr_imtime tau_mesh{params_.beta, Boson, params_.dlr_wmax, params_.dlr_eps, true};
    results->chiAB_tau = gf<mesh::dlr_imtime, tensor_valued<1>>{tau_mesh, make_shape(params.chi_ops.size())};
    chiAB_tau_.rebind(results->chiAB_tau.value());
    chiAB_tau_() = 0;

    // Precompute tau points from the DLR mesh
    L_ = tau_mesh.size();
    tau_points_.reserve(L_);
    for (auto tau : tau_mesh) tau_points_.push_back(tau_t::from_double(double(tau)));

    groups_ = finalize_groups(group_map, L_, n_blocks_);
  }

  void chiAB_tau::accumulate(mc_weight_t sign) {
    Z += sign;

    for (auto &grp : groups_) {
      long E = static_cast<long>(grp.terms.size());

      // Fill scratch arrays by calling each term with (tau_A, tau_B)
      for (long l = 0; l < L_; ++l) {
        std::array<tau_t, 2> taus = {tau_points_[l], tau_zero_};
        for (long e = 0; e < E; ++e) {
          auto const &blocks = grp.terms[e](std::span<const tau_t>{taus});
          for (int b = 0; b < n_blocks_; ++b) {
            int k = grp.signature[b];
            for (int p = 0; p < k; ++p) {
              grp.scratches[b].cs(l, e, p)    = blocks[b].cs[p];
              grp.scratches[b].cdags(l, e, p) = blocks[b].cdags[p];
            }
          }
        }
      }

      auto total_ratio = compute_insertion_ratios(grp, L_, qmc_config.dets, n_blocks_);

      // Scatter into chiAB_tau_
      for (long l = 0; l < L_; ++l)
        for (long e = 0; e < E; ++e) chiAB_tau_.data()(l, grp.terms[e].target_idx) += grp.terms[e].coef * sign * total_ratio(l, e);
    }
  }

  void chiAB_tau::collect_results(mpi::communicator const &comm) {
    Z          = mpi::all_reduce(Z, comm);
    chiAB_tau_ = mpi::all_reduce(chiAB_tau_, comm);
    chiAB_tau_ = chiAB_tau_ / Z;
  }

} // namespace triqs_ctint::measures
