// Copyright (c) 2018--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./histogram.hpp"

namespace triqs_ctint::measures {

  histogram::histogram(params_t const &, qmc_config_t const &qmc_config_, container_set *results)
     : qmc_config(qmc_config_), histogram_(results->histogram) {
    results->histogram = std::vector<double>(4);
  }

  void histogram::accumulate(mc_weight_t) {
    int k = qmc_config.perturbation_order();
    while (k >= histogram_->size()) histogram_->resize(2 * histogram_->size());
    histogram_.value()[k] += 1.;
    ++N;
  }

  void histogram::collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);

    // Make sure that all mpi threads have an equally sized histogram
    auto max_size = mpi::all_reduce(histogram_->size(), comm, MPI_MAX);
    histogram_->resize(max_size, 0.0);

    // Reduce histogram over all mpi threads
    histogram_ = mpi::all_reduce(histogram_.value(), comm);
    for (auto &h_k : histogram_.value()) h_k = h_k / N;
  }
} // namespace triqs_ctint::measures
