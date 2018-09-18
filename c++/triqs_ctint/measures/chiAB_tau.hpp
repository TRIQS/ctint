#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $\chi_{AB}(\tau)$ by operator insertion
  */
  struct chiAB_tau {

    chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results);

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters FIXME We cannot choose const, as we call try_insert
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t &qmc_config;

    // Container for the accumulation
    gf_view<imtime> chiAB_tau_;

    // The bosonic operator vectors
    using op_term_t = std::tuple<dcomplex, std::pair<int, int>, std::pair<cdag_t, c_t>>;
    std::vector<std::vector<op_term_t>> A_vec;
    std::vector<std::vector<op_term_t>> B_vec;

    // The average sign
    mc_weight_t Z = 0.0;

    // The tau-mesh
    gf_mesh<imtime> tau_mesh;
  };

} // namespace triqs_ctint::measures
