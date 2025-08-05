#pragma once
#include "../qmc_config.hpp"
#include "../container_set.hpp"
#include "triqs_ctint/types.hpp"
#include <triqs/stat/accumulator.hpp>

namespace triqs_ctint::measures {

  /**
  * Measure of $M_{ab}(\tau)$
  *
  * $M$ is the "reducible self-energy", see Eq. (41) in the Implementation Notes
  */
  struct M_tau_cheb {

    M_tau_cheb(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    /// Convert all accumulated samples to coefficients
    void convert_samples_to_coeffs();

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    std::vector<nda::matrix<std::vector<double>>> &tau_samples;
    std::vector<nda::matrix<std::vector<dcomplex>>> &weight_samples;
    std::vector<nda::array<dcomplex, 3>> &curlyG;
    std::vector<nda::array<double, 3>> auto_corr_times;
    std::vector<nda::array<triqs::stat::accumulator<dcomplex>, 3>> curlyG_acc;

    // Matrix views for the hartree term accumulation
    std::vector<matrix_view<M_tau_scalar_t>> M_hartree_;

    // The average sign
    mc_weight_t Z = 0.0;
  };

} // namespace triqs_ctint::measures
