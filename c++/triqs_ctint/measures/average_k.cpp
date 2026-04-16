// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./average_k.hpp"

namespace triqs_ctint::measures {

  average_k::average_k(params_t const &, qmc_config_t const &qmc_config_, container_set *results)
     : qmc_config(qmc_config_), average_k_(results->average_k), average_k_error_(results->average_k_error), k_bins_(dcomplex{0.0}, 128, 1) {
    average_k_ = 0.0;
  }

  void average_k::accumulate(mc_weight_t) {
    auto k = qmc_config.perturbation_order();
    average_k_ += k;
    k_bins_ << dcomplex(double(k));
    ++N;
  }

  void average_k::collect_results(mpi::communicator const &comm) {
    average_k_ = mpi::all_reduce(average_k_, comm);
    N          = mpi::all_reduce(N, comm);
    average_k_ = average_k_ / N;

    auto [m, err, tau] = k_bins_.mean_error_and_tau(comm);
    average_k_error_   = std::abs(err);
  }

  std::string average_k::report() const {
    std::ostringstream os;
    os << "Average perturbation order: " << average_k_ / N;
    return os.str();
  }

} // namespace triqs_ctint::measures
