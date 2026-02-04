// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include "../nfft_buf.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /**
  * Measure of $M^3_{abcd}(i\omega_1, i\omega_2)$
  *
  * $M^3$ is the essential building block for the fermion-boson vertices
  */
  struct M3ph_iw {

    M3ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_);

    // M3ph_iw needs to be uncopyable due to nfft_buf_t
    M3ph_iw(M3ph_iw const &)            = delete;
    M3ph_iw(M3ph_iw &&)                 = default;
    ~M3ph_iw()                          = default;
    M3ph_iw &operator=(M3ph_iw const &) = delete;
    M3ph_iw &operator=(M3ph_iw &&)      = delete;

    /// Accumulate M3ph measurement
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    chi3_dlr2d_iw_v_t M3ph_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers. buf_arrarr(block)(u_i,u_j)
    array<array<nfft_buf_t<2>, 2>, 1> buf_arrarr;
    array<array<nfft_buf_t<1>, 2>, 1> buf_arrarr_GM;
    array<array<nfft_buf_t<1>, 2>, 1> buf_arrarr_MG;

    // The non-interacting Green function
    g_tau_cv_t G0_tau;

    // Intermediate scattering matrix M on DLR2D mesh (type3 NFFT)
    block_gf<dlr2d_imfreq, matrix_valued> M;
    std::vector<std::array<mesh::matsubara_freq, 2>> target_mf_M;

    // Intermediate scattering matrices GM, MG on uniform imfreq mesh (type1 NFFT)
    block_gf<imfreq, matrix_valued> GM;
    block_gf<imfreq, matrix_valued> MG;
    array<array<dcomplex, 2>, 1> GMG;
  };

} // namespace triqs_ctint::measures
