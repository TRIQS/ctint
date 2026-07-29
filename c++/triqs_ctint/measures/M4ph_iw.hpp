// Copyright (c) 2023--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../qmc_config.hpp"
#include <triqs/utility/nfft/buffer.hpp>
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  // Same accumulation as iw4_accumulate, over the ph mesh (iW, iw, iwp):
  //
  //   M4(i,j,k,l) += sign * M1a(j,i) * M2a(l,k) - sign * M1b(l,i) * M2b(j,k)   [second term bl1 == bl2 only]
  //
  // Defined in M4ph_iw.cpp, declared here so benchmarks/iw_accum/iw4.cpp can drive it per block size.
  void iw4ph_accumulate(mc_weight_t sign, M4_M_t const &M, chi4_iw_v_t &M4ph, int bl1, int bl2, long bl2_size) noexcept;

  /**
  * Measure of $M^4_{abcd}(\tau_a, \tau_b, \tau_c)$
  *
  * $M^4$ is the essential building block for the two-particle Green function
  */
  struct M4ph_iw {

    M4ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results);

    // M4_iw needs to be uncopyable due to nfft::buffer_t
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
    chi4_iw_v_t M4ph_iw_;

    // The average sign
    mc_weight_t Z = 0.0;

    // Container of nfft_buffers: buf_arrarr(block)(u_j,u_i)
    array<array<nfft::buffer_t<2>, 2>, 1> buf_arrarr;

    // Intermediate scattering matrix in the measurement of M4
    M4_M_t M;
  };

} // namespace triqs_ctint::measures
