#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {
  struct f4ph_loc_iw {
    f4ph_loc_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_iw_t const &GinvG0_iw_, g_iw_t const &G0Ginv_iw_);

    // f4ph_loc_iw needs to be uncopyable due to nfft_buf_t
    f4ph_loc_iw(f4ph_loc_iw const &)            = delete;
    f4ph_loc_iw(f4ph_loc_iw &&)                 = default;
    ~f4ph_loc_iw()                              = default;
    f4ph_loc_iw &operator=(f4ph_loc_iw const &) = delete;
    f4ph_loc_iw &operator=(f4ph_loc_iw &&)      = delete;

    /// Accumulate M_tau using binning
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // Capture Ginv * G0 / G0 * Ginv
    g_iw_t const &GinvG0_iw;
    g_iw_t const &G0Ginv_iw;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    block2_gf_view<prod<imfreq, imfreq, imfreq>, tensor_valued<1>> f4ph_loc_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers: buf_arrarr(block)(u_j,u_i)
    array<array<nfft_buf_t<2>, 2>, 1> buf_arrarr;

    // Intermediate scattering matrix in the measurement of M4
    block_gf<prod<imfreq, imfreq>, matrix_valued> M;
    block_gf<prod<imfreq, imfreq>, tensor_valued<1>> m; // buffer for Ginv * G0 * M * G0 * Ginv
  };

} // namespace triqs_ctint::measures
