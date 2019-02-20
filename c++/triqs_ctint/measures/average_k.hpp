#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /// Measure of the average perturbation order
  struct average_k {

    average_k(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate average sign
    void accumulate(mc_weight_t sign);

    /// Reduce and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Reference to double for accumulation
    double &average_k_;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Accumulation counter
    long long N = 0;
  };

} // namespace triqs_ctint::measures
