#include "./vertex.hpp"

namespace triqs_ctint {

  double tau_t::beta = 0.0;

  // Note: The subtration of unsigned integers is inherently cyclic
  double cyclic_difference(tau_t const &tau1, tau_t const &tau2) { return double(tau_t{tau1.n - tau2.n}); }

} // namespace triqs_ctint
