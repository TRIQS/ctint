#include "./vertex.hpp"

namespace triqs_ctint {

  double tau_t::beta = 0.0;

  // Note: The subtration of unsigned integers is inherently cyclic
  double cyclic_difference(tau_t const &tau1, tau_t const &tau2) { return double(tau_t{tau1.n - tau2.n}); }

  tau_t make_tau_t(double tau) {
#ifdef DEBUG_CTINT
    if (tau < 0.0 or tau_t::beta < tau) TRIQS_RUNTIME_ERROR << " Tau-value outside [0,beta) interval not allowed in make_tau_t\n";
#endif
    return tau_t{uint32_t(tau_t::n_max / tau_t::beta * tau)};
  }

} // namespace triqs_ctint
