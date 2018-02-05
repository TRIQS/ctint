#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  // Measure the histogram of perturbation order
  struct histogram {

    histogram(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate perturbation order into histogram
    void accumulate(mc_weight_t sign);

    /// Reduce and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Reference to accumulation vector
    std::vector<double> &histogram_;

    // Accumulation counter
    int N = 0;
  };

} // namespace triqs_ctint::measures
