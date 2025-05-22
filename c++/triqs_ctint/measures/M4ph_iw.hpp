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
  struct M4ph_iw {

    M4ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    // M4_iw needs to be uncopyable due to nfft_buf_t
    M4ph_iw(M4ph_iw const &)            = delete;
    M4ph_iw(M4ph_iw &&)                 = default;
    ~M4ph_iw()                          = default;
    M4ph_iw &operator=(M4ph_iw const &) = delete;
    M4ph_iw &operator=(M4ph_iw &&)      = delete;

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
    block2_gf_view<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> M4ph_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers: buf_arrarr(block)(u_j,u_i)
    array<array<nfft_buf_t<2>, 2>, 1> buf_arrarr;

    // Intermediate scattering matrix in the measurement of M4
    using M_layout = nda::contiguous_layout_with_stride_order<nda::encode(std::array{0, 1, 3, 2})>;
    block_gf<prod<imfreq, imfreq>, matrix_valued, M_layout> M;
  };

} // namespace triqs_ctint::measures
