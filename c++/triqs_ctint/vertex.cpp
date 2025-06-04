// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./vertex.hpp"

namespace triqs_ctint {

  double tau_t::beta = 0.0;

  // Note: The subtration of unsigned integers is inherently cyclic
  std::pair<double, double> cyclic_difference(tau_t const &tau1, tau_t const &tau2) {
    // Assume tau1 > tau2 for the equal time case
    double sign  = tau2 > tau1 ? -1.0 : 1.0;
    double value = double(tau_t{tau1.n - tau2.n});
    return std::make_pair(sign, value);
  }

  std::pair<double, double> cyclic_difference(double tau1, double tau2) {
    // Assume tau1 > tau2 for the equal time case
    double dtau  = tau1 - tau2;
    int nshifts  = static_cast<int>(std::floor(dtau / tau_t::beta));
    double sign  = (nshifts % 2 == 0) ? 1.0 : -1.0;
    double value = dtau - nshifts * tau_t::beta;
    return std::make_pair(sign, value);
  }

  tau_t make_tau_t(double tau) {
#ifdef DEBUG_CTINT
    if (tau < 0.0 or tau_t::beta < tau) TRIQS_RUNTIME_ERROR << " Tau-value outside [0,beta) interval not allowed in make_tau_t\n";
#endif
    return tau_t{uint32_t(tau_t::n_max / tau_t::beta * tau)};
  }

  std::ostream &operator<<(std::ostream &os, vertex_idx_t const &v) {
    os << "vertex_idx_t{"
       << "b1,u1=" << v.b1 << "," << v.u1 << ", "
       << "b2,u2=" << v.b2 << "," << v.u2 << ", "
       << "b3,u3=" << v.b3 << "," << v.u2 << ", "
       << "b4,u4=" << v.b4 << "," << v.u2 << "}";
    return os;
  }

  std::ostream &operator<<(std::ostream &os, vertex_t const &v) {
    os << "vertex_t{"
       << "idx=" << v.idx << ", "
       << "tau1=" << v.tau1 << ", "
       << "tau2=" << v.tau2 << ", "
       << "tau3=" << v.tau3 << ", "
       << "tau4=" << v.tau4 << ", "
       << "amplitude=" << v.amplitude << ", "
       << "proposition_proba=" << v.proposition_proba << ", "
       << "s=" << v.s << "}";
    return os;
  }

} // namespace triqs_ctint
