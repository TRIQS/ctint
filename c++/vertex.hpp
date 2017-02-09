#pragma once
#include "./types.hpp"
#include <triqs/det_manip/det_manip.hpp>

namespace triqs_ctint {

  /**
   * A point in imaginary time, i.e. $\tau \in [0,\beta]$, but defined on a very fine grid.
   * The position in the segment is given by an uint64_t, i.e. a very long integer. 
   * This allows exact comparisons, which notoriously dangerous on floating point number.
   */
  struct tau_t {

    /// Maximum value that can be stored inside a uint64_t
    static constexpr uint64_t n_max = std::numeric_limits<uint64_t>::max();

    /// Inverse temperature associated with all $\tau$ points
    static double beta;

    /// $\tau$ value, represented as an integer on a very fine grid
    uint64_t n = 0;

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
  };

  /// Function returns the result of the difference of two tau points, shifted to the interval [0,\beta]
  double cyclic_difference(tau_t const &tau1, tau_t const &tau2);

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
    double amplitude;

    /// Probability of proposition for this vertex
    double proposition_proba;

    /// Switch for alpha shift
    bool has_alpha_shift = false;

    /// Value of auxiliary spin
    int s = 0;
  };

} // namespace triqs_ctint
