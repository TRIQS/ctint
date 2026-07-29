// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include <triqs/utility/nfft/buffer.hpp>
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  // The pp update over the sparse DLR2D mesh:
  //
  //   M3(i,j,k,l) += sign * M1a(j,i) * M2a(l,k) - sign * M1b(l,i) * M2b(j,k)   [second term bl1 == bl2 only]
  //
  // Defined in M3pp_iw.cpp.
  void dlr2d_iw3pp_accumulate(mc_weight_t sign, M3_G_t const &GM, chi3_dlr2d_iw_v_t M3, int bl1, int bl2, long bl2_size) noexcept;

  /**
  * Measure of $M^3_{abcd}(i\omega_1, i\omega_2)$
  *
  * $M^3$ is the essential building block for the fermion-boson vertices
  */
  struct M3pp_iw {

    M3pp_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_);

    // M3pp_iw needs to be uncopyable due to nfft::buffer_t
    M3pp_iw(M3pp_iw const &)            = delete;
    M3pp_iw(M3pp_iw &&)                 = default;
    ~M3pp_iw()                          = default;
    M3pp_iw &operator=(M3pp_iw const &) = delete;
    M3pp_iw &operator=(M3pp_iw &&)      = delete;

    /// Accumulate M3pp measurement
    void accumulate(mc_weight_t sign);

    /// Collect results and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Capture the parameters
    params_t const &params;

    // The Monte-Carlo configuration
    qmc_config_t const &qmc_config;

    // Container for the accumulation
    chi3_dlr2d_iw_v_t M3pp_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers. buf_arrarr(block)(u_i,u_j)
    array<array<nfft::buffer_t<1>, 2>, 1> buf_arrarr;

    // The non-interacting Green function
    g_tau_cv_t G0_tau;

    // Intermediate scattering matrix using type1 NFFT on uniform imfreq mesh
    M3_G_t GM;
  };

} // namespace triqs_ctint::measures
