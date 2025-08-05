#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"
#include "triqs_ctint/types.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M_{ab}(\tau)$
  *
  * $M$ is the "reducible self-energy", see Eq. (41) in the Implementation Notes
  */
  struct M_iw_dlr {

    M_iw_dlr(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    void convert_samples_to_dlr();

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    std::vector<nda::matrix<std::vector<double>>> &tau_samples;
    std::vector<nda::matrix<std::vector<dcomplex>>> &weight_samples;
    block_gf<dlr_imfreq, matrix_valued> &M_iw_dlr_;
    // Matrix views for the hartree term accumulation
    std::vector<matrix_view<M_tau_scalar_t>> M_hartree_;

    // The average sign
    mc_weight_t Z = 0.0;
  };

} // namespace triqs_ctint::measures
