#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^\sigma_{ab}(\tau)$
  *
  * $M$ is the "reducible self-energy", see Eq. (31) in the Implementation Notes
  */
  struct M_tau {

    M_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(double sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    block_gf_view<imtime, matrix_valued> M_tau_;

    // The average sign
    double Z = 0.0;
  };

} // namespace triqs_ctint::measures
