// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include <optional>
#include <span>

namespace triqs_ctint::measures {

  // Per-block operator arguments with tau set by call operator
  struct block_ops_t {
    std::vector<c_t> cs;
    std::vector<cdag_t> cdags;
  };

  // Per-operator configuration
  struct op_config_t {
    int orbital;
    tau_t tau_offset; // added to the base tau selected by tau_idx
    int tau_idx;    // selects base tau from the vector passed to operator()
  };

  // Per-block configuration
  struct block_config_t {
    std::vector<op_config_t> cdag_config;
    std::vector<op_config_t> c_config;
  };

  // Block signature: k per block (0 = uninvolved). Length n_blocks.
  using block_signature_t = std::vector<int>;

  // A single decomposed monomial term with callable tau-filling
  struct operator_term_t {
    long target_idx; // index into the output array:
                     //   static_obs: which operator in the static_obs list
                     //   chiAB_tau:  which (A,B) pair in chi_ops
    dcomplex coef;   // includes sign from Wick contraction and block factorization

    // Produce per-block c_t/cdag_t with tau set.
    // Caller passes ALL base taus: each op selects one via its tau_idx.
    std::vector<block_ops_t> const &operator()(std::span<const tau_t> taus) const {
      for (size_t b = 0; b < config_.size(); ++b) {
        auto const &cfg = config_[b];
        auto &blk       = blocks_[b];
        for (size_t p = 0; p < cfg.c_config.size(); ++p) {
          auto const &c = cfg.c_config[p];
          blk.cs[p]     = c_t{taus[c.tau_idx] + c.tau_offset, c.orbital};
        }
        for (size_t p = 0; p < cfg.cdag_config.size(); ++p) {
          auto const &cd = cfg.cdag_config[p];
          blk.cdags[p]   = cdag_t{taus[cd.tau_idx] + cd.tau_offset, cd.orbital};
        }
      }
      return blocks_;
    }

    // Convenience: single tau (covers static_obs)
    std::vector<block_ops_t> const &operator()(tau_t tau) const { return operator()(std::span<const tau_t>{&tau, 1}); }

    block_signature_t signature() const {
      block_signature_t sig(config_.size());
      for (size_t b = 0; b < config_.size(); ++b) sig[b] = static_cast<int>(config_[b].c_config.size());
      return sig;
    }

    // Construction helper: only used by make_static_term / make_chiAB_term
    static operator_term_t make(long target_idx, dcomplex coef, std::vector<block_config_t> config, int n_blocks) {
      operator_term_t t;
      t.target_idx = target_idx;
      t.coef       = coef;
      t.config_    = std::move(config);
      t.blocks_.resize(n_blocks);
      for (int b = 0; b < n_blocks; ++b) {
        t.blocks_[b].cs.resize(t.config_[b].c_config.size());
        t.blocks_[b].cdags.resize(t.config_[b].cdag_config.size());
      }
      return t;
    }

    private:
    mutable std::vector<block_ops_t> blocks_; // pre-allocated output, length n_blocks
    std::vector<block_config_t> config_;      // length n_blocks
  };

  // Per-block scratch arrays for batched det_manip calls
  struct block_scratch_t {
    nda::array<c_t, 3> cs;       // (L, E, k)
    nda::array<cdag_t, 3> cdags; // (L, E, k)
  };

  // Group of terms sharing the same block signature
  struct term_group_t {
    block_signature_t signature;
    std::vector<operator_term_t> terms;
    std::vector<block_scratch_t> scratches; // length n_blocks
  };

  // ==================== Sign computation ====================

  // Count inversions in a sequence (O(n^2), n is typically <= 6)
  inline int count_inversions(std::vector<int> const &seq) {
    int inv = 0;
    for (size_t i = 0; i < seq.size(); ++i)
      for (size_t j = i + 1; j < seq.size(); ++j)
        if (seq[i] > seq[j]) ++inv;
    return inv;
  }

  // Compute the sign relating the product of per-block determinants to the physical
  // expectation value of a normally-ordered operator.
  //
  // Given a canonical monomial c†_1...c†_n c_n...c_1 with tau offsets assigned as
  // c's at 0..n-1 (lower), c†'s at n..2n-1 (higher), the Wick theorem gives:
  //
  //   <op> = (-1)^{n(n-1)/2 + P_row + P_col} × product_of_block_dets
  //
  // where P_row = inversions to sort c's by block, P_col = inversions to sort c†'s by block.
  //
  // cdag_blocks: block index of each c†, in canonical monomial order
  // c_blocks: block index of each c, in canonical monomial order
  inline int compute_factorization_sign(std::vector<int> const &cdag_blocks, std::vector<int> const &c_blocks) {
    int n     = static_cast<int>(cdag_blocks.size());
    int P_col = count_inversions(cdag_blocks);
    int P_row = count_inversions(c_blocks);
    int total = n * (n - 1) / 2 + P_row + P_col;
    return (total % 2 == 0) ? 1 : -1;
  }

  // ==================== Factory: static_obs term ====================

  // Decompose a monomial into an operator_term_t where all ops use a single tau (tau_idx = 0).
  // Tau offsets: c's at 0..n-1 (lower), c†'s at n..2n-1 (higher).
  // Returns nullopt if any block has unbalanced #cdag != #c.
  inline std::optional<std::pair<operator_term_t, block_signature_t>> make_static_term(monomial_t const &m, dcomplex coef, long target_idx,
                                                                                       hilbert_space::gf_struct_t const &gf_struct, int n_blocks) {
    using triqs::operators::get_int_indices;

    if (m.size() % 2 != 0) return std::nullopt; // odd degree
    int n = static_cast<int>(m.size()) / 2;    // number of c†-c pairs

    // Extract block indices for c†'s and c's in canonical monomial order.
    // Canonical: c†'s come first (all dagger), then c's (all non-dagger).
    std::vector<int> cdag_blocks, c_blocks;
    std::vector<int> cdag_orbitals, c_orbitals;
    for (auto const &op : m) {
      auto [bl, orb] = get_int_indices(op, gf_struct);
      if (op.dagger) {
        cdag_blocks.push_back(bl);
        cdag_orbitals.push_back(orb);
      } else {
        c_blocks.push_back(bl);
        c_orbitals.push_back(orb);
      }
    }
    if (static_cast<int>(cdag_blocks.size()) != n || static_cast<int>(c_blocks.size()) != n) return std::nullopt;

    // Check per-block balance
    std::vector<int> n_cdag_per_block(n_blocks, 0), n_c_per_block(n_blocks, 0);
    for (int bl : cdag_blocks) ++n_cdag_per_block[bl];
    for (int bl : c_blocks) ++n_c_per_block[bl];
    for (int b = 0; b < n_blocks; ++b)
      if (n_cdag_per_block[b] != n_c_per_block[b]) return std::nullopt;

    // Build per-block configs with tau offsets.
    // c's get LOWER offsets 0..n-1, c†'s get HIGHER offsets n..2n-1.
    // This ensures cdag.tau > c.tau, which makes G0hat evaluate at tau ≈ β
    // giving f(c, cdag) = <c† c> (the density) for equal-time operators.
    std::vector<block_config_t> config(n_blocks);
    for (int i = 0; i < n; ++i) config[c_blocks[i]].c_config.push_back({c_orbitals[i], tau_t{std::uint64_t(i)}, 0});
    for (int i = 0; i < n; ++i) config[cdag_blocks[i]].cdag_config.push_back({cdag_orbitals[i], tau_t{std::uint64_t(n + i)}, 0});

    // Compute sign
    int sign = compute_factorization_sign(cdag_blocks, c_blocks);

    auto term = operator_term_t::make(target_idx, coef * static_cast<double>(sign), std::move(config), n_blocks);
    auto sig  = term.signature();
    return std::make_pair(std::move(term), std::move(sig));
  }

  // ==================== Factory: chiAB term ====================

  // Decompose a combined A*B monomial into an operator_term_t.
  // A-ops use tau_idx = 0 (variable tau), B-ops use tau_idx = 1 (caller passes tau_t::get_zero()).
  // Tau offsets: c's at 0..n-1 (lower), c†'s at n..2n-1 (higher), with A/B interleaved within each group.
  // Returns nullopt if any block has unbalanced #cdag != #c.
  inline std::optional<std::pair<operator_term_t, block_signature_t>> make_chiAB_term(monomial_t const &A_mono, monomial_t const &B_mono,
                                                                                      dcomplex coef, long target_idx,
                                                                                      hilbert_space::gf_struct_t const &gf_struct, int n_blocks) {
    using triqs::operators::get_int_indices;

    // Collect all operators from A and B, separated into c†'s and c's.
    // Sort c†'s ascending by (block, orbital) to match TRIQS canonical normal ordering.
    // Sort c's descending by (block, orbital) to match TRIQS canonical normal ordering.
    struct op_info_t {
      int block, orbital, tau_idx;
    };
    std::vector<op_info_t> cdags, cs;

    for (auto const &op : B_mono) {
      auto [bl, orb] = get_int_indices(op, gf_struct);
      if (op.dagger)
        cdags.push_back({bl, orb, 1});
      else
        cs.push_back({bl, orb, 1});
    }
    for (auto const &op : A_mono) {
      auto [bl, orb] = get_int_indices(op, gf_struct);
      if (op.dagger)
        cdags.push_back({bl, orb, 0});
      else
        cs.push_back({bl, orb, 0});
    }

    // Sort c†'s ascending (matching canonical c†-ordering in TRIQS)
    std::stable_sort(cdags.begin(), cdags.end(),
                     [](auto const &a, auto const &b) { return std::tie(a.block, a.orbital) < std::tie(b.block, b.orbital); });
    // Sort c's descending (matching canonical c-ordering in TRIQS: c's in reverse index order)
    std::stable_sort(cs.begin(), cs.end(), [](auto const &a, auto const &b) { return std::tie(a.block, a.orbital) > std::tie(b.block, b.orbital); });

    int n = static_cast<int>(cdags.size());
    if (n != static_cast<int>(cs.size())) return std::nullopt;

    // Check per-block balance
    std::vector<int> n_cdag_per_block(n_blocks, 0), n_c_per_block(n_blocks, 0);
    for (auto const &op : cdags) ++n_cdag_per_block[op.block];
    for (auto const &op : cs) ++n_c_per_block[op.block];
    for (int b = 0; b < n_blocks; ++b)
      if (n_cdag_per_block[b] != n_c_per_block[b]) return std::nullopt;

    // Build per-block configs. c's at LOWER offsets 0..n-1, c†'s at HIGHER offsets n..2n-1.
    std::vector<block_config_t> config(n_blocks);
    for (int i = 0; i < n; ++i) config[cs[i].block].c_config.push_back({cs[i].orbital, tau_t{std::uint64_t(i)}, cs[i].tau_idx});
    for (int i = 0; i < n; ++i) config[cdags[i].block].cdag_config.push_back({cdags[i].orbital, tau_t{std::uint64_t(n + i)}, cdags[i].tau_idx});

    // === Analytical sign: (-1)^{n_B_cdag * m_A_c + P_row + P_col} ===
    //
    // Physical: <A(tau) B(0)> with A's operators interleaving B's.
    // Moving B's n_B_cdag c-daggers past A's m_A_c c-operators: (-1)^{n_B_cdag * m_A_c}.
    // Then Wick theorem with specific row/col orderings.
    //
    // Wick c† order: [A's c†'s ascending, B's c†'s ascending] by (block, orbital)
    // Wick c order: [B's c's ascending, A's c's ascending] by (block, orbital)
    // Config c† order: block-grouped, ascending within block (= cdags list order)
    // Config c order: block-grouped, descending within block (= cs list, distributed to blocks)
    //
    // P_row = inversions from Wick c†-order to config c†-order
    // P_col = inversions from Wick c-order to config c-order

    // Count A's c-daggers and c-operators separately (they differ for unbalanced monomials)
    auto wick_cdags = cdags;
    auto mid        = std::stable_partition(wick_cdags.begin(), wick_cdags.end(), [](auto const &op) { return op.tau_idx == 0; });
    int n_B_cdag    = n - static_cast<int>(std::distance(wick_cdags.begin(), mid));
    int m_A_c       = static_cast<int>(std::count_if(cs.begin(), cs.end(), [](auto const &op) { return op.tau_idx == 0; }));

    // Build Wick c ordering: [B's ascending, A's ascending]
    auto asc     = [](auto const &a, auto const &b) { return std::tie(a.block, a.orbital) < std::tie(b.block, b.orbital); };
    auto wick_cs = cs;
    std::sort(wick_cs.begin(), wick_cs.end(), asc);
    std::stable_partition(wick_cs.begin(), wick_cs.end(), [](auto const &op) { return op.tau_idx == 1; });

    // Config c† order = cdags (already block-grouped ascending).
    // Config c order = cs grouped by block, descending within each block.
    auto config_cs = cs;
    std::stable_sort(config_cs.begin(), config_cs.end(), [](auto const &a, auto const &b) { return a.block < b.block; });

    // P_row/P_col: inversions from Wick ordering to config ordering.
    // Build permutation by matching (block, orbital, tau_idx), then count inversions.
    auto build_permutation = [](auto const &from, auto const &to) {
      int sz = static_cast<int>(from.size());
      std::vector<bool> used(sz, false);
      std::vector<int> perm(sz);
      for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
          if (!used[j] && from[i].block == to[j].block && from[i].orbital == to[j].orbital && from[i].tau_idx == to[j].tau_idx) {
            perm[i] = j;
            used[j] = true;
            break;
          }
        }
      }
      return perm;
    };

    auto perm_row = build_permutation(wick_cdags, cdags);
    auto perm_col = build_permutation(wick_cs, config_cs);

    int P_row = count_inversions(perm_row);
    int P_col = count_inversions(perm_col);
    int total = n_B_cdag * m_A_c + P_row + P_col;
    int sign  = (total % 2 == 0) ? 1 : -1;

    auto term = operator_term_t::make(target_idx, coef * static_cast<double>(sign), std::move(config), n_blocks);
    auto sig  = term.signature();
    return std::make_pair(std::move(term), std::move(sig));
  }

  // ==================== Shared utilities ====================

  // Flatten group_map into a vector and pre-allocate scratch arrays for each group.
  inline std::vector<term_group_t> finalize_groups(std::map<block_signature_t, term_group_t> &group_map, long L, int n_blocks) {
    std::vector<term_group_t> groups;
    groups.reserve(group_map.size());
    for (auto &[sig, grp] : group_map) {
      long E = static_cast<long>(grp.terms.size());
      grp.scratches.resize(n_blocks);
      for (int b = 0; b < n_blocks; ++b) {
        int k = sig[b];
        if (k > 0) {
          grp.scratches[b].cs.resize(L, E, k);
          grp.scratches[b].cdags.resize(L, E, k);
        }
      }
      groups.push_back(std::move(grp));
    }
    return groups;
  }

  // Compute the (L, E) array of total insertion ratios for one term group.
  inline nda::array<dcomplex, 2> compute_insertion_ratios(term_group_t &grp, long L, std::vector<det_t> &dets, int n_blocks) {
    long E = static_cast<long>(grp.terms.size());
    nda::array<dcomplex, 2> total_ratio(L, E);
    total_ratio() = 1.0;

    for (int b = 0; b < n_blocks; ++b) {
      int k = grp.signature[b];
      if (k == 0) continue;
      auto &det = dets[b];
      auto &sc  = grp.scratches[b];

      if (k == 1) {
        total_ratio *= det.insert_ratios(0, 0, sc.cs(nda::range::all, nda::range::all, 0), sc.cdags(nda::range::all, nda::range::all, 0));
      } else if (k == 2) {
        total_ratio *= det.insert2_ratios(0, 1, 0, 1, sc.cs(nda::range::all, nda::range::all, 0), sc.cs(nda::range::all, nda::range::all, 1),
                                          sc.cdags(nda::range::all, nda::range::all, 0), sc.cdags(nda::range::all, nda::range::all, 1));
      } else {
        for (long l = 0; l < L; ++l)
          total_ratio(l, nda::range::all) *= det.insertk_ratios(sc.cs(l, nda::range::all, nda::range::all), sc.cdags(l, nda::range::all, nda::range::all));
      }
    }
    return total_ratio;
  }

} // namespace triqs_ctint::measures
