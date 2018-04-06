#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
   * Measure of $M^\sigma_{ab}(i\omega)$
   *
   * $M$ is the "reducible self-energy", see Eq. (31) in the Implementation Notes
   */
  struct M_iw {

    M_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    // M_iw needs to be uncopyable due to nfft_buf_t
    M_iw(M_iw const &) = delete;
    M_iw(M_iw &&)      = default;
    ~M_iw()            = default;
    M_iw &operator=(M_iw const &) = delete;
    M_iw &operator=(M_iw &&) = default;

    /// Accumulate M_iw using nfft
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    block_gf_view<imfreq, matrix_valued> M_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers. buf_vec[block_idx](a,b)
    std::vector<array<nfft_buf_t<1>, 2>> buf_vec;
  };

} // namespace triqs_ctint::measures
