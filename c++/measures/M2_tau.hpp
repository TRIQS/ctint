#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^\2_{abcd}(\tau)$
  *
  * $M^2$ is the essential building block for the susceptibilities
  */
  template <Chan_t Chan> struct M2_tau {

    M2_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, block_gf<imtime, matrix_valued> const &G0_tau_);

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
    block2_gf_view<imtime, tensor_valued<4>> M2_tau_;

    // The average sign
    double Z = 0.0;

    // The non-interacting Green function
    block_gf<imtime, matrix_valued> const &G0_tau;
  };

} // namespace triqs_ctint::measures

#include "M2_tau_impl.hpp"
