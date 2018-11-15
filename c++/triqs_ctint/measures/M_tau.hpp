#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M_{ab}(\tau)$
  *
  * $M$ is the "reducible self-energy", see Eq. (41) in the Implementation Notes
  */
  struct M_tau {

    M_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Gf view for the M_tau accumulation
    block_gf_view<imtime, M_tau_target_t> M_tau_;

    // Matrix views for the hartree term accumulation
    std::vector<matrix_view<M_tau_scalar_t>> M_hartree_;

    // The average sign
    mc_weight_t Z = 0.0;
  };

} // namespace triqs_ctint::measures
