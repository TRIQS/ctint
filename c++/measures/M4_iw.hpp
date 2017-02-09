#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^\4_{abcd}(\tau_a, \tau_b, \tau_c)$
  *
  * $M^4$ is the essential building block for the two-particle Green function
  */
  struct M4_iw {

    using block_type = gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>;
    using view_type  = block2_gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>;

    M4_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    // M4_iw needs to be uncopyable due to nfft_buf_t
    M4_iw(M4_iw const &) = delete;
    M4_iw(M4_iw &&)      = default;
    M4_iw &operator=(M4_iw const &) = delete;
    M4_iw &operator=(M4_iw &&) = default;

    /// Accumulate M_tau using binning
    void accumulate(double sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    view_type M4_iw_;

    // The average sign
    double Z = 0.0;

    // Container of nfft_buffers. buf_arr(block_1,block_2)(a,b)
    std::vector<std::vector<array<triqs::utility::nfft_buf_t<3>, 4>>> buf_vecvec; // FIXME New block_array type?
  };

} // namespace triqs_ctint::measures
