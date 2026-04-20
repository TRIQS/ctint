// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include <triqs/utility/nfft/buffer.hpp>
#include <triqs/utility/nfft/matrix_buffer.hpp>
#include "../container_set.hpp"


namespace triqs_ctint::measures {

  /**
  * Measure of $M^3_{abcd}(i\omega_1, i\omega_2)$
  *
  * $M^3$ is the essential building block for the fermion-boson vertices
  */
  struct M3ph_iw {

    M3ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_);

    // M3ph_iw needs to be uncopyable due to nfft::buffer_t
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

    // Container of nfft_buffers for GM and MG (type1 NFFT, uniform grid)
    array<array<nfft::buffer_t<1>, 2>, 1> buf_arrarr_GM;
    array<array<nfft::buffer_t<1>, 2>, 1> buf_arrarr_MG;

    // The non-interacting Green function
    g_tau_cv_t G0_tau;

    // Intermediate scattering matrix M on DLR2D mesh
    using M_layout = nda::contiguous_layout_with_stride_order<nda::encode(std::array{0, 2, 1})>;
    block_gf<dlr2d_imfreq, matrix_valued, M_layout> M;

    // Matrix NFFT buffer for M (factored product-grid DFT), one per block
    std::vector<nfft::matrix_buffer_t<M_layout>> M_bufs;

    // Intermediate scattering matrices GM, MG on uniform imfreq mesh (type1 NFFT)
    using GM_layout = nda::contiguous_layout_with_stride_order<nda::encode(std::array{0, 2, 1})>;
    block_gf<imfreq, matrix_valued, GM_layout> GM;
    using MG_layout = nda::contiguous_layout_with_stride_order<nda::encode(std::array{0, 2, 1})>;
    block_gf<imfreq, matrix_valued, MG_layout> MG;
    array<array<dcomplex, 2, nda::F_layout>, 1> GMG;

  };

} // namespace triqs_ctint::measures
