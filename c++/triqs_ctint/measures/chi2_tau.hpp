#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $\chi^2_{abcd}(\tau)$ by operator insertion
  */
  template <Chan_t Chan> struct chi2_tau {

    chi2_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters FIXME We cannot choose const, as we call try_insert
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t &qmc_config;

    // Container for the accumulation
    block2_gf_view<imtime, tensor_valued<4>> chi2_tau_;

    // The average sign
    mc_weight_t Z = 0.0;

    // The tau-mesh
    mesh::imtime tau_mesh;
  };

} // namespace triqs_ctint::measures

#include "./chi2_tau_impl.hpp"
