// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./average_k.hpp"

namespace triqs_ctint::measures {

  average_k::average_k(params_t const &, qmc_config_t const &qmc_config_, container_set *results)
     : qmc_config(qmc_config_), average_k_(results->average_k) {
    average_k_ = 0.0;
  }

  void average_k::accumulate(mc_weight_t) {
    average_k_ += qmc_config.perturbation_order();
    ++N;
  }

  void average_k::collect_results(mpi::communicator const &comm) {
    average_k_ = mpi::all_reduce(average_k_, comm);
    N          = mpi::all_reduce(N, comm);
    average_k_ = average_k_ / N;
  }

  std::string average_k::report() const {
    std::ostringstream os;
    os << "Average perturbation order: " << average_k_ / N;
    return os.str();
  }

} // namespace triqs_ctint::measures
