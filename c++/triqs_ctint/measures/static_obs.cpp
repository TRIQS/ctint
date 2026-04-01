// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./static_obs.hpp"

namespace triqs_ctint::measures {

  static_obs::static_obs(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_),
       qmc_config(qmc_config_),
       results_(results),
       n_blocks_(static_cast<int>(params_.gf_struct.size())),
       L_(params_.n_tau_static_obs) {

    auto const &ops = params.static_obs;
    if (ops.empty()) TRIQS_RUNTIME_ERROR << "Empty operator list in static_obs measurement";
    long n_obs = ops.size();

    // Init result accumulator and per-step contribution array
    result_.resize(n_obs);
    result_() = 0;
    step_contrib_.resize(n_obs);
    step_contrib_() = 0;

    // Init linear binning accumulators for error analysis (one per observable)
    obs_bins_.reserve(n_obs);
    for (long i = 0; i < n_obs; ++i) obs_bins_.emplace_back(dcomplex{0.0}, 128, 1);

    // Init constant parts
    constant_parts_ = nda::array<dcomplex, 1>(n_obs);
    constant_parts_() = 0;

    if (L_ <= 0) TRIQS_RUNTIME_ERROR << "n_tau_static_obs must be positive, got " << L_;

    // Precompute uniform tau grid: tau_k = k * beta / L for k = 0, ..., L-1
    tau_points_.resize(L_);
    for (long k = 0; k < L_; ++k) tau_points_[k] = tau_t::from_double(k * params.beta / L_);

    // Unified decomposition: all monomials go through make_static_term
    std::map<block_signature_t, term_group_t> group_map;

    for (auto [obs_idx, C] : itertools::enumerate(ops)) {
      for (auto const &term : C) {
        auto const &m = term.monomial;

        if (m.empty()) {
          constant_parts_(obs_idx) += dcomplex(term.coef);
        } else if (m.size() % 2 != 0) {
          TRIQS_RUNTIME_ERROR << "Odd-degree monomial (degree " << m.size() << ") in static_obs";
        } else {
          auto result = make_static_term(m, dcomplex(term.coef), static_cast<long>(obs_idx), params.gf_struct, n_blocks_);
          if (!result) continue; // block-unbalanced -> zero contribution
          auto &[oterm, sig]       = *result;
          group_map[sig].signature = sig;
          group_map[sig].terms.push_back(std::move(oterm));
        }
      }
    }

    groups_ = finalize_groups(group_map, L_, n_blocks_);
  }

  void static_obs::accumulate(mc_weight_t sign) {
    Z += sign;
    ++N_;
    step_contrib_() = 0;

    for (auto &grp : groups_) {
      long E = static_cast<long>(grp.terms.size());

      // Fill scratch arrays by calling each term
      for (long l = 0; l < L_; ++l) {
        auto tau = tau_points_[l];
        for (long e = 0; e < E; ++e) {
          auto const &blocks = grp.terms[e](tau);
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

      // Scatter into step_contrib_
      for (long e = 0; e < E; ++e) {
        auto val = grp.terms[e].coef * sign;
        for (long l = 0; l < L_; ++l) step_contrib_(grp.terms[e].target_idx) += val * total_ratio(l, e);
      }
    }

    // Constant parts: no tau dependence, multiply by L_ to compensate the /L_ normalization in collect_results
    for (long i = 0; i < constant_parts_.size(); ++i) step_contrib_(i) += constant_parts_(i) * sign * L_;

    // Accumulate into result and feed per-step contributions into bins
    result_ += step_contrib_;
    for (long i = 0; i < step_contrib_.size(); ++i) obs_bins_[i] << step_contrib_(i);
  }

  void static_obs::collect_results(mpi::communicator const &comm) {
    Z       = mpi::all_reduce(Z, comm);
    N_      = mpi::all_reduce(N_, comm);
    result_ = mpi::all_reduce(result_, comm);
    result_ = result_ / (Z * L_);
    results_->static_obs = result_;

    // Compute error bars from linear binning
    nda::array<double, 1> errors(result_.size());
    for (long i = 0; i < result_.size(); ++i) {
      auto [m, err, tau] = obs_bins_[i].mean_error_and_tau(comm);
      errors(i)          = std::abs(err) / (std::abs(Z / N_) * L_);
    }
    results_->static_obs_errors = errors;
  }

} // namespace triqs_ctint::measures
