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
  struct M3pp_tau {

    M3pp_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Gf view for the accumulation of M3pp_tau
    block2_gf_view<prod<imtime, imtime>, tensor_valued<4>> M3pp_tau_;

    // Gf view for the accumulation of M3pp_delta
    block2_gf_view<imtime, tensor_valued<4>> M3pp_delta_;

    // The average sign
    mc_weight_t Z = 0.0;

    // The non-interacting Green function
    g_tau_cv_t G0_tau;

    // The tau-mesh for a single argument
    mesh::imtime tau_mesh;

    // A helper type to precompute the tau-binning
    struct idx_t {
      int tau_idx;
      int u;
      tau_t tau_pt;
    };
  };

} // namespace triqs_ctint::measures
