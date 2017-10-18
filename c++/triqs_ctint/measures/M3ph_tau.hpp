#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^3_{abcd}(\tau_1, \tau_2)$
  *
  * $M^3$ is the essential building block for the fermion-boson verticies
  */
  struct M3ph_tau {

    M3ph_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, block_gf<imtime, matrix_valued> const &G0_tau_);

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
    block2_gf_view<cartesian_product<imtime, imtime>, tensor_valued<4>> M3ph_tau_;

    // The average sign
    double Z = 0.0;

    // The non-interacting Green function
    block_gf<imtime, matrix_valued> const &G0_tau;

    // The tau-mesh for a single argument
    gf_mesh<imtime> tau_mesh;

    // A helper type to precompute the tau-binning
    struct idx_t {
      int tau_idx;
      int u;
    };
  };

} // namespace triqs_ctint::measures
