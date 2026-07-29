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

  // The ph update over the sparse DLR2D mesh:
  //
  //   M3(i,j,k,l) += sign * M1a(j,i) * M2a(l,k) - sign * M1b(l,i) * M2b(j,k)   [second term bl1 == bl2 only]
  //
  // Defined in M3ph_iw.cpp.
  void dlr2d_iw3ph_accumulate(mc_weight_t sign, M3_M_t const &M, M3_GMG_t const &GMG, M3_G_t const &GM, M3_G_t const &MG, chi3_dlr2d_iw_v_t &M3,
                              int bl1, int bl2, long bl2_size) noexcept;

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

    // Scattering matrix M on DLR2D mesh (factored product-grid DFT)
    M3_M_t M;
    std::vector<nfft::matrix_buffer_t<M3_M_layout>> M_bufs;

    // Intermediate scattering matrices GM, MG on uniform imfreq mesh (type1 NFFT)
    M3_G_t GM;
    M3_G_t MG;
    M3_GMG_t GMG;
  };

} // namespace triqs_ctint::measures
