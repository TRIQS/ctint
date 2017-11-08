// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "./types.hpp"
#include <triqs/det_manip/det_manip.hpp>

namespace triqs_ctint {

  /**
   * A point in imaginary time, i.e. $\tau \in [0,\beta]$, but defined on a very fine grid.
   * The position in the segment is given by an uint32_t, i.e. a very long integer.
   * This allows exact comparisons, which notoriously dangerous on floating point number.
   */
  struct tau_t {

    /// Maximum value that can be stored inside a uint32_t
    static constexpr uint32_t n_max = std::numeric_limits<uint32_t>::max();

    /// Inverse temperature associated with all $\tau$ points
    static double beta;

    /// $\tau$ value, represented as an integer on a very fine grid
    uint32_t n = 0;

    /// Get a random point in $[0,\beta[$
    template <typename RNG> static tau_t get_random(RNG &rng) { return tau_t{rng(n_max)}; }

    /// Cast to corresponding double value in $[0,\beta]$
    explicit operator double() const { return beta * double(n) / n_max; }

    // --- Comparison operators
    bool operator==(const tau_t &tau) const { return n == tau.n; }
    bool operator!=(const tau_t &tau) const { return n != tau.n; }
    bool operator<(const tau_t &tau) const { return n < tau.n; }
    bool operator>(const tau_t &tau) const { return n > tau.n; }
    bool operator<=(const tau_t &tau) const { return n <= tau.n; }
    bool operator>=(const tau_t &tau) const { return n >= tau.n; }

    /// Operator allowing output to std::ostream
    friend std::ostream &operator<<(std::ostream &out, tau_t const &tau) {
      return out << double(tau) << " [tau_t : beta = " << tau.beta << " n = " << tau.n << "]";
    }

    /// Return \tau = 0
    static constexpr tau_t get_zero() { return tau_t{0}; }

    /// Return \tau = 0^{+} = 0 + \delta
    static constexpr tau_t get_zero_plus() { return tau_t{1}; }

    /// Return \tau = 0^{++} = 0 + 2*\delta
    static constexpr tau_t get_zero_plus_plus() { return tau_t{2}; }

    /// Return \tau = \beta
    static constexpr tau_t get_beta() { return tau_t{n_max}; }

    /// Return \tau = beta^{-} = \beta - \delta
    static constexpr tau_t get_beta_minus() { return tau_t{n_max - 1}; }

    /// Return \tau = beta^{--} = \beta - 2*\delta
    static constexpr tau_t get_beta_minus_minus() { return tau_t{n_max - 2}; }
  };

  /// Calculate the time-difference of two tau points shifted to the interval [0,\beta] as well
  /// as the sign change resulting from the shift in a fermionic function
  std::pair<double, double> cyclic_difference(tau_t const &tau1, tau_t const &tau2);
  std::pair<double, double> cyclic_difference(double tau1, double tau2);

  // Generate a tau_t-object from a double
  tau_t make_tau_t(double tau);

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
