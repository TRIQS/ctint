// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "./types.hpp"
#include <triqs/utility/tau_t.hpp>
#include <triqs/det_manip/det_manip.hpp>

namespace triqs_ctint {

  /// Imaginary-time type, provided by TRIQS.
  using tau_t = triqs::utility::tau_t;

  /// Bring cyclic_difference into this namespace.
  using triqs::utility::cyclic_difference;

  /**
   * Type representing the set of discrete quantum numbers for the vertices of
   * the microscopic model at hand. We distinguish between block-indeces
   * (in which the bare Green function is diagonal) and non-block indeces.
   */
  struct vertex_idx_t {

    /// First operator of the vertex (outgoing, c^\dagger): block index, non-block
    int b1, u1;

    /// Second operator of the vertex (ingoing, c): block index, non-block
    int b2, u2;

    /// Third operator of the vertex (outgoing, c^\dagger): block index, non-block
    int b3, u3;

    /// Fourth operator of the vertex (ingoing, c): block index, non-block
    int b4, u4;
  };

  std::ostream &operator<<(std::ostream &os, vertex_idx_t const &v);

  /**
   * Type representing an interaction vertex of the microscopic model at hand.
   * Can be inserted in the Monte-Carlo move.
   */
  struct vertex_t {

    /// Object containing discrete quantum numbers for external legs, i.e. block and non-block index
    vertex_idx_t idx;

    /// Imaginary times for all four external legs (c^\dagger, c, c^\dagger, c)
    tau_t tau1, tau2, tau3, tau4;

    /// Amplitude of the vertex, i.e. U, U(tau1-tau2), etc...
    U_scalar_t amplitude;

    /// Probability of proposition for this vertex
    double proposition_proba;

    /// The label of the vertex (position in h_int)
    int vertex_label = 0;

    /// Value of auxiliary spin
    int s = 0;
  };

  std::ostream &operator<<(std::ostream &os, vertex_t const &v);

} // namespace triqs_ctint
