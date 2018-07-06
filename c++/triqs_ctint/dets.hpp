#pragma once

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

    /// Switch for alpha shift along the diagonal of the matrix. Relevant for density-type interactions
    bool with_alpha_shift = false;

    /// The auxiliary spin
    int s = 0;

    /// Lexicographical sorting of arg_t. This determines the order of row and columns inside the dets.
    bool operator<(arg_t const &x) const { return std::tie(tau, u, s) < std::tie(x.tau, x.u, x.s); }
  };

  using c_t    = arg_t<false>;
  using cdag_t = arg_t<true>;

  /**
   * Functor that evaluates the matrix elements, used by the det_manip. 
   * Cares for the possible alpha-shift along the diagonal of the matrix.
   */
  struct G0hat_t {
    /// The (shifted) non-interacting Green function
    gf_const_view<imtime, g_tau_t::target_t> G0_shift_tau;

    /// The alpha function
    array<double, 2> alpha;

    g_tau_t::target_t::scalar_t operator()(c_t const &c, cdag_t const &cdag) const {
      if ((c.tau == cdag.tau) and (c.u == cdag.u)) { return G0_shift_tau[0](c.u, cdag.u) + (c.with_alpha_shift ? 1.0 - alpha(c.u, c.s) : 1.0); }
      double d_tau = cyclic_difference(c.tau, cdag.tau);
      auto res     = G0_shift_tau[closest_mesh_pt(d_tau)](c.u, cdag.u);
      // For the equal-time case, consider the order <c cdag>
      return (c.tau >= cdag.tau ? res : -res);
    }
  };

  /// Type of a single determinant
  using det_t = det_manip<G0hat_t>;

} // namespace triqs_ctint
