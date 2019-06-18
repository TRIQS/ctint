#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^4_{abcd}(\tau_a, \tau_b, \tau_c)$
  *
  * $M^4$ is the essential building block for the two-particle Green function
  */
  struct M4_iw {

    M4_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    // M4_iw needs to be uncopyable due to nfft_buf_t
    M4_iw(M4_iw const &) = delete;
    M4_iw(M4_iw &&)      = default;
    ~M4_iw()             = default;
    M4_iw &operator=(M4_iw const &) = delete;
    M4_iw &operator=(M4_iw &&) = delete;

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    block2_gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> M4_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers: buf_arrarr(block)(u_j,u_i)
    array<array<nfft_buf_t<2>, 2>, 1> buf_arrarr;

    // Intermediate scattering matrix in the measurement of M4
    block_gf<cartesian_product<imfreq, imfreq>, matrix_valued> M;
  };

} // namespace triqs_ctint::measures
