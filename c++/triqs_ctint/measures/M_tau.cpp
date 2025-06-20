// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./M_tau.hpp"

namespace triqs_ctint::measures {

  M_tau::M_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    // Init measurement container and capture view
    results->M_tau = block_gf<imtime, M_tau_target_t>{{params.beta, Fermion, params.n_tau}, params.gf_struct};
    M_tau_.rebind(results->M_tau.value());
    M_tau_() = 0;

    results->M_hartree = make_block_vector<M_tau_scalar_t>(params.gf_struct);
    for (auto &m : results->M_hartree.value()) M_hartree_.push_back(m);
  }

  void M_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    // for (auto const & [D,M] : triqs::std::zip(config.dets, results->M_tau)) 	// C++17
    for (int b = 0; b < M_tau_.size(); ++b) {
      // Loop over every index pair (x,y) in the determinant matrix
      // for (auto const & [x,y,Ginv] : D ) 	// C++17
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {
        // Check for the equal-time case
        if (c_i.tau == cdag_j.tau) {
          M_hartree_[b](cdag_j.u, c_i.u) += Ginv * sign;
        } else {

          // Absolut time-difference tau of the index pair
          auto [s, dtau] = cyclic_difference(cdag_j.tau, c_i.tau);

          // Project tau to closest point on the binning grid
          M_tau_[b][closest_mesh_pt(dtau)](cdag_j.u, c_i.u) += Ginv * s * sign;
        }
      });
    }
  }

  void M_tau::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z           = mpi::all_reduce(Z, comm);
    M_tau_      = mpi::all_reduce(M_tau_, comm);
    double dtau = M_tau_[0].mesh().delta();
    M_tau_      = M_tau_ / (-Z * dtau * params.beta);
    for (auto &m : M_hartree_) {
      m = mpi::all_reduce(m, comm);
      m = m / (-Z * params.beta);
    }

    // Multiply first and last bin by two
    for (auto &M : M_tau_) {
      M[0] *= 2.0;
      M[params.n_tau - 1] *= 2.0;
    }
  }

} // namespace triqs_ctint::measures
