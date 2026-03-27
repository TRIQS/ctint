// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./static_obs.hpp"
#include "./../types.hpp"
#include <triqs/operators/util/extractors.hpp>

using namespace triqs::utility;

namespace triqs_ctint::measures {

  static_obs::static_obs(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), results_(results), L_(params_.n_tau_static_obs) {

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

    // Maps for grouping
    std::map<int, bilinear_group_t> bilinear_map;
    using quartic_key_t = std::tuple<quartic_case_t, int, int>;
    std::map<quartic_key_t, quartic_group_t> quartic_map;

    for (auto [obs_idx, C] : itertools::enumerate(ops)) {
      for (auto const &term : C) {
        auto const &m = term.monomial;

        if (m.size() == 0) {
          // Constant term
          constant_parts_(obs_idx) += dcomplex(term.coef);

        } else if (m.size() == 2) {
          // Bilinear c†c
          if (!m[0].dagger or m[1].dagger)
            TRIQS_RUNTIME_ERROR << "Bilinear monomial in static_obs not of the form c^+ c";
          auto [bl_cdag, idx_cdag] = get_int_indices(m[0], params.gf_struct);
          auto [bl_c, idx_c]       = get_int_indices(m[1], params.gf_struct);
          if (bl_cdag != bl_c)
            TRIQS_RUNTIME_ERROR << "Off-diagonal block bilinear in static_obs not supported (bl_cdag=" << bl_cdag << ", bl_c=" << bl_c << ")";
          auto &grp  = bilinear_map[bl_cdag];
          grp.bl_det = bl_cdag;
          grp.entries.push_back({static_cast<long>(obs_idx), dcomplex(term.coef), idx_cdag, idx_c});

        } else if (m.size() == 4) {
          // Quartic c†c†cc (normal ordered)
          if (!(m[0].dagger && m[1].dagger && !m[2].dagger && !m[3].dagger))
            TRIQS_RUNTIME_ERROR << "Quartic monomial in static_obs not of the form c^+ c^+ c c";

          // Pair A: (m[0], m[3]) = (c†_a, c_b), Pair B: (m[1], m[2]) = (c†_c, c_d)
          auto [bl_cdag_A, idx_cdag_A] = get_int_indices(m[0], params.gf_struct);
          auto [bl_c_A, idx_c_A]       = get_int_indices(m[3], params.gf_struct);
          auto [bl_cdag_B, idx_cdag_B] = get_int_indices(m[1], params.gf_struct);
          auto [bl_c_B, idx_c_B]       = get_int_indices(m[2], params.gf_struct);

          auto info = classify_quartic_blocks(bl_cdag_A, bl_c_A, bl_cdag_B, bl_c_B, dcomplex(term.coef));
          if (!info) continue;

          quartic_key_t key{info->case_type, info->bl_det1, info->bl_det2};
          auto &grp     = quartic_map[key];
          grp.case_type = info->case_type;
          grp.bl_det1   = info->bl_det1;
          grp.bl_det2   = info->bl_det2;
          grp.entries.push_back({static_cast<long>(obs_idx), info->coef, idx_cdag_A, idx_c_A, idx_cdag_B, idx_c_B});

        } else {
          TRIQS_RUNTIME_ERROR << "Operator degree " << m.size() << " not supported in static_obs (max 4)";
        }
      }
    }

    // Flatten maps to vectors and pre-allocate scratch arrays
    bilinear_groups_.reserve(bilinear_map.size());
    for (auto &[key, grp] : bilinear_map) {
      long E = static_cast<long>(grp.entries.size());
      grp.cs.resize(L_, E);
      grp.cdags.resize(L_, E);
      bilinear_groups_.push_back(std::move(grp));
    }
    quartic_groups_.reserve(quartic_map.size());
    for (auto &[key, grp] : quartic_map) {
      long E = static_cast<long>(grp.entries.size());
      grp.c_A.resize(L_, E);
      grp.c_B.resize(L_, E);
      grp.cdag_A.resize(L_, E);
      grp.cdag_B.resize(L_, E);
      quartic_groups_.push_back(std::move(grp));
    }
  }

  void static_obs::accumulate(mc_weight_t sign) {
    Z += sign;
    ++N_;
    step_contrib_() = 0;

    // Scatter ratios into per-observable step contributions (e-outer for entry lookup hoisting)
    auto accum = [&](auto const &entries, auto const &ratios) {
      long E = static_cast<long>(entries.size());
      for (long e = 0; e < E; ++e) {
        auto const &entry = entries[e];
        auto val          = entry.coef * sign;
        for (long l = 0; l < L_; ++l) step_contrib_(entry.target_idx) += val * ratios(l, e);
      }
    };

    // --- Bilinear groups ---
    for (auto &grp : bilinear_groups_) {
      auto &det = qmc_config.dets[grp.bl_det];
      for (long l = 0; l < L_; ++l) {
        auto tau      = tau_points_[l];
        auto tau_plus = tau + tau_t::epsilon();
        for (long e = 0; e < static_cast<long>(grp.entries.size()); ++e) {
          auto const &entry = grp.entries[e];
          grp.cs(l, e)     = c_t{tau, entry.idx_c};
          grp.cdags(l, e)  = cdag_t{tau_plus, entry.idx_cdag};
        }
      }
      accum(grp.entries, det.insert_ratios(0, 0, grp.cs, grp.cdags));
    }

    // --- Quartic groups ---
    for (auto &grp : quartic_groups_) {
      // Fill pre-allocated (L_, E) arrays for all four operators
      // Distinct tau offsets (+0,+1,+2,+3) enforce strict time-ordering for the det_manip insertions
      for (long l = 0; l < L_; ++l) {
        auto tau = tau_points_[l];
        for (long e = 0; e < static_cast<long>(grp.entries.size()); ++e) {
          auto const &entry  = grp.entries[e];
          grp.c_B(l, e)     = c_t{tau, entry.idx_c_B};
          grp.cdag_B(l, e)  = cdag_t{tau + tau_t{std::uint64_t{1}}, entry.idx_cdag_B};
          grp.c_A(l, e)     = c_t{tau + tau_t{std::uint64_t{2}}, entry.idx_c_A};
          grp.cdag_A(l, e)  = cdag_t{tau + tau_t{std::uint64_t{3}}, entry.idx_cdag_A};
        }
      }

      switch (grp.case_type) {
        case quartic_case_t::AAAA: {
          accum(grp.entries, qmc_config.dets[grp.bl_det1].insert2_ratios(0, 1, 0, 1, grp.c_A, grp.c_B, grp.cdag_A, grp.cdag_B));
          break;
        }
        case quartic_case_t::AABB: {
          auto r1 = qmc_config.dets[grp.bl_det1].insert_ratios(0, 0, grp.c_A, grp.cdag_A);
          auto r2 = qmc_config.dets[grp.bl_det2].insert_ratios(0, 0, grp.c_B, grp.cdag_B);
          r1 *= r2;
          accum(grp.entries, r1);
          break;
        }
        case quartic_case_t::ABAB: {
          auto r1 = qmc_config.dets[grp.bl_det1].insert_ratios(0, 0, grp.c_B, grp.cdag_A);
          auto r2 = qmc_config.dets[grp.bl_det2].insert_ratios(0, 0, grp.c_A, grp.cdag_B);
          r1 *= r2;
          accum(grp.entries, r1);
          break;
        }
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
