#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $G^2_{abcc}(\tau,0,0,0)$ by operator insertion
  */
  struct G2_fluct_diag_tau {

    G2_fluct_diag_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    /// Accumulate G2_fluct_diag_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters FIXME We cannot choose const, as we call try_insert
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t &qmc_config;

    // Container for the accumulation
    block2_gf_view<imtime, tensor_valued<4>> G2_fluct_diag_tau_;

    // The average sign
    mc_weight_t Z = 0.0;

    // The tau-mesh
    gf_mesh<imtime> tau_mesh;
  };

} // namespace triqs_ctint::measures

#include "./G2_fluct_diag_tau_impl.hpp"
