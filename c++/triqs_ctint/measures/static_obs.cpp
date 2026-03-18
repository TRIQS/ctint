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
    for (long k = 0; k < L_; ++k) tau_points_[k] = make_tau_t(k * params.beta / L_);

    // Maps for grouping
    std::map<int, bilinear_group_t> bilinear_map;                                       // keyed by bl_det
    using quartic_key_t = std::tuple<case_t, int, int>;
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

          bool is_AABB = (bl_cdag_A == bl_c_A && bl_cdag_B == bl_c_B);
          bool is_ABAB = (bl_cdag_A == bl_c_B && bl_cdag_B == bl_c_A);
          bool is_AAAA = is_AABB && is_ABAB;

          case_t ct;
          int bl_det1, bl_det2;
          dcomplex coef = dcomplex(term.coef);

          if (is_AAAA) {
            ct = case_t::AAAA;
            bl_det1 = bl_cdag_A;
            bl_det2 = bl_cdag_A;
          } else if (is_AABB) {
            ct = case_t::AABB;
            bl_det1 = bl_cdag_A;
            bl_det2 = bl_cdag_B;
          } else if (is_ABAB) {
            ct = case_t::ABAB;
            bl_det1 = bl_cdag_A;
            bl_det2 = bl_c_A;
            coef = -coef; // minus sign from operator swap
          } else {
            continue; // other block combinations vanish
          }

          quartic_key_t key{ct, bl_det1, bl_det2};
          auto &grp     = quartic_map[key];
          grp.case_type = ct;
          grp.bl_det1   = bl_det1;
          grp.bl_det2   = bl_det2;
          grp.entries.push_back({static_cast<long>(obs_idx), coef, idx_cdag_A, idx_c_A, idx_cdag_B, idx_c_B});

        } else {
          TRIQS_RUNTIME_ERROR << "Operator degree " << m.size() << " not supported in static_obs (max 4)";
        }
      }
    }

    // Flatten maps to vectors
    bilinear_groups_.reserve(bilinear_map.size());
    for (auto &[key, grp] : bilinear_map) bilinear_groups_.push_back(std::move(grp));
    quartic_groups_.reserve(quartic_map.size());
    for (auto &[key, grp] : quartic_map) quartic_groups_.push_back(std::move(grp));
  }

  void static_obs::accumulate(mc_weight_t sign) {
    Z += sign;
    ++N_;
    step_contrib_() = 0;

    // --- Bilinear groups ---
    for (auto const &grp : bilinear_groups_) {
      long E     = static_cast<long>(grp.entries.size());
      auto &det  = qmc_config.dets[grp.bl_det];

      nda::array<c_t, 2> cs(L_, E);
      nda::array<cdag_t, 2> cdags(L_, E);

      for (long l = 0; l < L_; ++l) {
        auto tau = tau_points_[l];
        auto tau_plus = tau_t{tau.n + 1};
        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          cs(l, e)    = c_t{tau, entry.idx_c};
          cdags(l, e) = cdag_t{tau_plus, entry.idx_cdag};
        }
      }

      auto ratios = det.insert_ratios(0, 0, cs, cdags); // shape (L_, E)

      for (long l = 0; l < L_; ++l)
        for (long e = 0; e < E; ++e) step_contrib_(grp.entries[e].obs_idx) += grp.entries[e].coef * sign * ratios(l, e);
    }

    // --- Quartic groups ---
    for (auto const &grp : quartic_groups_) {
      long E = static_cast<long>(grp.entries.size());

      // Fill rank-2 (L_, E) arrays for all four operators — common to all cases
      nda::array<c_t, 2> c_A(L_, E), c_B(L_, E);
      nda::array<cdag_t, 2> cdag_A(L_, E), cdag_B(L_, E);
      for (long l = 0; l < L_; ++l) {
        auto tau = tau_points_[l];
        for (long e = 0; e < E; ++e) {
          auto const &entry = grp.entries[e];
          c_B(l, e)    = c_t{tau_t{tau.n}, entry.idx_c_B};
          cdag_B(l, e) = cdag_t{tau_t{tau.n + 1}, entry.idx_cdag_B};
          c_A(l, e)    = c_t{tau_t{tau.n + 2}, entry.idx_c_A};
          cdag_A(l, e) = cdag_t{tau_t{tau.n + 3}, entry.idx_cdag_A};
        }
      }

      // Accumulate ratios into step_contrib_
      auto accum = [&](auto const &ratios) {
        for (long l = 0; l < L_; ++l)
          for (long e = 0; e < E; ++e) step_contrib_(grp.entries[e].obs_idx) += grp.entries[e].coef * sign * ratios(l, e);
      };

      // Compute ratios — only the insert call differs between cases
      switch (grp.case_type) {
        case case_t::AAAA: {
          accum(qmc_config.dets[grp.bl_det1].insert2_ratios(0, 1, 0, 1, c_A, c_B, cdag_A, cdag_B));
          break;
        }
        case case_t::AABB: {
          auto r1 = qmc_config.dets[grp.bl_det1].insert_ratios(0, 0, c_A, cdag_A);
          auto r2 = qmc_config.dets[grp.bl_det2].insert_ratios(0, 0, c_B, cdag_B);
          accum(r1 * r2);
          break;
        }
        case case_t::ABAB: {
          auto r1 = qmc_config.dets[grp.bl_det1].insert_ratios(0, 0, c_B, cdag_A);
          auto r2 = qmc_config.dets[grp.bl_det2].insert_ratios(0, 0, c_A, cdag_B);
          accum(r1 * r2);
          break;
        }
      }
    }

    // --- Constant parts (no tau dependence, multiply by L_ to compensate collect_results /L_) ---
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
