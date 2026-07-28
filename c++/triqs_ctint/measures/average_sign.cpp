// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./average_sign.hpp"

namespace triqs_ctint::measures {

  average_sign::average_sign(params_t const &, qmc_config_t const &, container_set *results)
     : average_sign_(results->average_sign),
       average_sign_error_(results->average_sign_error),
       nmeasures(results->nmeasures),
       sign_bins_(dcomplex{0.0}, 128, 1) {
    average_sign_ = 0.0;
    nmeasures     = 0;
  }

  void average_sign::accumulate(mc_weight_t sign) {
    average_sign_ += sign;
    sign_bins_ << dcomplex(sign);
    ++count;
  }

  void average_sign::collect_results(mpi::communicator const &comm) {
    average_sign_ = mpi::all_reduce(average_sign_, comm);
    count         = mpi::all_reduce(count, comm);
    average_sign_ = average_sign_ / count;
    nmeasures     = count;

    auto [m, err, tau]  = sign_bins_.mean_error_and_tau(comm);
    average_sign_error_ = std::abs(err);
  }

  std::string average_sign::report() const {
    std::ostringstream os;
    os << "Average sign: " << average_sign_ / count;
    return os.str();
  }

} // namespace triqs_ctint::measures
