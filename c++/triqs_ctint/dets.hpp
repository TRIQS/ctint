// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include <ostream>

namespace triqs_ctint {

  using triqs::det_manip::det_manip;

  //------------------------------------

  /**
   * Type of row and column argument of the Green function matrix inside the determinant.
   * G(x,y) is evaluated at the time-difference x.tau - y.tau. 
   */
  template <bool dag> struct arg_t {

    /// C (false) or Cdag (true)
    static constexpr bool dagger = dag;

    /// The imaginary time
    tau_t tau;

    /// The orbital (or non-block) index
    int u;

    /// The label of the vertex (position in h_int). Defaults to -1 for external operators.
    int vertex_label = -1;

    /// The position in the quartic operator (cdag_0 c_0 cdag_1 c_1). Defaults to -1 for external operators.
    int pos = -1;

    /// The auxiliary spin. Defaults to -1 for external operators.
    int s = -1;

    /// Lexicographical sorting of arg_t. This determines the order of row and columns inside the dets.
    bool operator<(arg_t const &x) const { return std::tie(tau, u, vertex_label, pos, s) < std::tie(x.tau, x.u, x.vertex_label, x.pos, x.s); }
  };

  using c_t    = arg_t<false>;
  using cdag_t = arg_t<true>;

  template <bool dag> std::ostream &operator<<(std::ostream &os, arg_t<dag> const &c) {
    os << (c.dagger ? "cdag_t{" : "c_t{") << "tau=" << c.tau << ", "
       << "u=" << c.u << ", "
       << "vertex_label=" << c.vertex_label << ", "
       << "pos=" << c.pos << ", "
       << "s=" << c.s << "}";
    return os;
  }

  /**
   * Functor that evaluates the matrix elements, used by the det_manip. 
   * Cares for the possible alpha-shift along the diagonal of the matrix.
   */
  struct G0hat_t {
    /// The (shifted) non-interacting Green function
    gf_const_view<imtime, g_tau_t::target_t> G0_shift_tau;

    /// The alpha function
    array_const_view<g_tau_scalar_t, 4> alpha;

    // Precomputed constants for fast mesh index computation
    double delta_inv_ = 0;
    long n_tau_       = 0;
    long n_orb_       = 0;

    G0hat_t(gf_const_view<imtime, g_tau_t::target_t> g0, array_const_view<g_tau_scalar_t, 4> a)
       : G0_shift_tau(std::move(g0)), alpha(std::move(a)) {
      auto const &mesh = G0_shift_tau.mesh();
      delta_inv_       = mesh.delta_inv();
      n_tau_           = mesh.size();
      n_orb_           = G0_shift_tau.data().shape()[1];
    }

    G0hat_t()                            = default;
    G0hat_t(G0hat_t const &)             = default;
    G0hat_t(G0hat_t &&)                  = default;
    G0hat_t &operator=(G0hat_t const &)  = delete;
    G0hat_t &operator=(G0hat_t &&)       = delete;

    g_tau_t::target_t::scalar_t operator()(c_t const &c, cdag_t const &cdag) const {
      // Equal-time contractions between operators of the same interaction vertex get an alpha shift.
      // External measurement operators carry vertex_label == -1 and must not receive an alpha shift.
      if (c.tau == cdag.tau) {
        TRIQS_ASSERT2(c.vertex_label == cdag.vertex_label,
                       "Equal-time operators must have the same vertex_label:\n  " << c << "\n  " << cdag);
        auto G0_val = G0_shift_tau.data()(0, c.u, cdag.u) + (c.u == cdag.u ? 1 : 0);
        return c.vertex_label < 0 ? G0_val : G0_val - alpha(c.vertex_label, cdag.pos, c.pos, c.s);
      }
      // Compute sign and dtau via cyclic_difference
      auto [sign, dtau] = cyclic_difference(c.tau, cdag.tau);
      // Compute mesh data index directly: idx = clamp(round(dtau * delta_inv), 0, n_tau-1)
      long idx     = std::clamp(static_cast<long>(dtau * delta_inv_ + 0.5), 0L, n_tau_ - 1);
      return sign * G0_shift_tau.data()(idx, c.u, cdag.u);
    }
  };

  /// Type of a single determinant
  using det_t = det_manip<G0hat_t>;

} // namespace triqs_ctint
