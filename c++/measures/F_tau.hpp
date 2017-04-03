#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  // Measure for the F = \Sigma G
  struct F_tau {

    F_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, block_gf<imtime, matrix_valued> const &G0_tau_);

    /// Accumulate F_tau
    void accumulate(double sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    block_gf_view<imtime, matrix_valued> F_tau_;

    // The average sign
    double Z = 0.0;

    // An access to the Green's functions
    block_gf<imtime, matrix_valued> const &G0_tau;
  };

} // namespace triqs_ctint::measures
