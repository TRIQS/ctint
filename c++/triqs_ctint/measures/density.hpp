#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of the density matrix by operator insertion
  */
  struct density {

    density(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters FIXME We cannot choose const, as we call try_insert
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t &qmc_config;

    // Container for the accumulation
    std::vector<matrix_view<dcomplex>> density_;

    // The average sign
    mc_weight_t Z = 0.0;
  };

} // namespace triqs_ctint::measures
