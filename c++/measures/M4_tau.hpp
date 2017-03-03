#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^4_{abcd}(\tau_a, \tau_b, \tau_c)$
  *
  * $M^4$ is the essential building block for the two-particle Green function
  */
  struct M4_tau {

    using block_type = gf<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>>;
    using view_type  = block2_gf_view<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>>;

    M4_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

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
    view_type M4_tau_;

    // The average sign
    double Z = 0.0;
  };

} // namespace triqs_ctint::measures
