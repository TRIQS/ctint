#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^3_{abcd}(i\omega_1, i\omega_2)$
  *
  * $M^3$ is the essential building block for the fermion-boson verticies
  */
  struct M3ph_iw {

    M3ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_);

    // M3ph_iw needs to be uncopyable due to nfft_buf_t
    M3ph_iw(M3ph_iw const &) = delete;
    M3ph_iw(M3ph_iw &&)      = default;
    ~M3ph_iw()               = default;
    M3ph_iw &operator=(M3ph_iw const &) = delete;
    M3ph_iw &operator=(M3ph_iw &&) = default;

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    block2_gf_view<cartesian_product<imfreq, imfreq>, tensor_valued<4>> M3ph_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers. buf_arrarr(block)(u_i,u_j)
    array<array<nfft_buf_t<2>, 2>, 1> buf_arrarr;
    array<array<nfft_buf_t<1>, 2>, 1> buf_arrarr_GM;
    array<array<nfft_buf_t<1>, 2>, 1> buf_arrarr_MG;

    // The non-interacting Green function
    g_tau_cv_t G0_tau;

    // Intermediate scattering matrix in the measurement of M3ph
    block_gf<cartesian_product<imfreq, imfreq>, matrix_valued> M;
    block_gf<imfreq, matrix_valued> GM;
    block_gf<imfreq, matrix_valued> MG;
    array<array<dcomplex, 2>, 1> GMG;
  };

} // namespace triqs_ctint::measures
