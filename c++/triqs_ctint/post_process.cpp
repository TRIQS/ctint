// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "post_process.hpp"
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

namespace triqs_ctint {

  chi4_iw_t G2_conn_from_M4(chi4_iw_t::const_view_type M4_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    chi4_iw_t G2_conn_iw = M4_iw; // FIXME Product Ranges with += Lazy Expressions

    double beta  = M_iw[0].mesh().beta();
    int n_blocks = M_iw.size();

    // Calculate connected part of M4
    chi4_iw_t M4_iw_conn = M4_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        M4_iw_conn(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << M4_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
              - beta * kronecker(iw1_, iw2_) * M_iw[bl1](iw1_)(j_, i_) * M_iw[bl2](iw3_)(l_, k_)
              + beta * kronecker(bl1, bl2) * kronecker(iw2_, iw3_) * M_iw[bl1](iw1_)(l_, i_) * M_iw[bl2](iw3_)(j_, k_);

    // Calculate disconnected part of the two-particle Green function
    G2_conn_iw() = 0.;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = M4_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M4_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                G2_conn_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << G2_conn_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
                      + G0_iw[bl1](iw2_)(j_, n) * G0_iw[bl2](iw1_ - iw2_ + iw3_)(l_, p) * M4_iw_conn(bl1, bl2)(iw1_, iw2_, iw3_)(m, n, o, p)
                         * G0_iw[bl1](iw1_)(m, i_) * G0_iw[bl2](iw3_)(o, k_);
      }

    return G2_conn_iw;
  }

  chi4_iw_t G2pp_conn_from_M4pp(chi4_iw_t::const_view_type M4pp_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    chi4_iw_t G2pp_conn_iw = M4pp_iw; // FIXME Product Ranges with += Lazy Expressions

    double beta  = M_iw[0].mesh().beta();
    int n_blocks = M_iw.size();

    // Calculate connected part of M4
    chi4_iw_t M4pp_iw_conn = M4pp_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        M4pp_iw_conn(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << M4pp_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
              - beta * kronecker(iw_, iW_ - iwp_) * M_iw[bl1](iw_)(j_, i_) * M_iw[bl2](iW_ - iw_)(l_, k_)
              + beta * kronecker(bl1, bl2) * kronecker(iW_ - iwp_, iW_ - iw_) * M_iw[bl1](iw_)(l_, i_) * M_iw[bl2](iW_ - iw_)(j_, k_);

    // Calculate disconnected part of the two-particle Green function
    G2pp_conn_iw() = 0.;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = M4pp_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M4pp_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                G2pp_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2pp_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + G0_iw[bl1](iW_ - iwp_)(j_, n) * G0_iw[bl2](iwp_)(l_, p) * M4pp_iw_conn(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p)
                         * G0_iw[bl1](iw_)(m, i_) * G0_iw[bl2](iW_ - iw_)(o, k_);
      }

    return G2pp_conn_iw;
  }

  chi4_iw_t G2ph_conn_from_M4ph(chi4_iw_t::const_view_type M4ph_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    chi4_iw_t G2ph_conn_iw = M4ph_iw; // FIXME Product Ranges with += Lazy Expressions

    double beta  = M_iw[0].mesh().beta();
    int n_blocks = M_iw.size();

    // Calculate connected part of M4
    chi4_iw_t M4ph_iw_conn = M4ph_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        M4ph_iw_conn(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << M4ph_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
              - beta * kronecker(iw_, iW_ + iw_) * M_iw[bl1](iw_)(j_, i_) * M_iw[bl2](iW_ + iwp_)(l_, k_)
              + beta * kronecker(bl1, bl2) * kronecker(iW_ + iw_, iW_ + iwp_) * M_iw[bl1](iw_)(l_, i_) * M_iw[bl2](iW_ + iwp_)(j_, k_);

    // Calculate disconnected part of the two-particle Green function
    G2ph_conn_iw() = 0.;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = M4ph_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M4ph_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                G2ph_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2ph_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + G0_iw[bl1](iW_ + iw_)(j_, n) * G0_iw[bl2](iwp_)(l_, p) * M4ph_iw_conn(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p)
                         * G0_iw[bl1](iw_)(m, i_) * G0_iw[bl2](iW_ + iwp_)(o, k_);
      }

    return G2ph_conn_iw;
  }

  chi4_iw_t F_from_G2c(chi4_iw_t::const_view_type G2_conn_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();

    // Temporary quantities
    g_iw_t Ginv = inverse(G_iw);

    // Calculate vertex function F
    chi4_iw_t F_iw = G2_conn_iw; // FIXME Product Ranges with += Lazy Expressions
    F_iw()         = 0;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = G2_conn_iw(bl1, bl2).target_shape()[0];
        int bl2_size = G2_conn_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                F_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << F_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
                      + Ginv[bl1](iw2_)(j_, n) * Ginv[bl2](iw1_ - iw2_ + iw3_)(l_, p) * G2_conn_iw(bl1, bl2)(iw1_, iw2_, iw3_)(m, n, o, p)
                         * Ginv[bl1](iw1_)(m, i_) * Ginv[bl2](iw3_)(o, k_);
      }

    return F_iw;
  }

  chi4_iw_t Fpp_from_G2pp_conn(chi4_iw_t::const_view_type G2pp_conn_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();

    // Temporary quantities
    g_iw_t Ginv = inverse(G_iw);

    // Calculate vertex function F
    chi4_iw_t Fpp_iw = G2pp_conn_iw; // FIXME Product Ranges with += Lazy Expressions
    Fpp_iw()         = 0;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = G2pp_conn_iw(bl1, bl2).target_shape()[0];
        int bl2_size = G2pp_conn_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                Fpp_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << Fpp_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + Ginv[bl1](iW_ - iwp_)(j_, n) * Ginv[bl2](iwp_)(l_, p) * G2pp_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p) * Ginv[bl1](iw_)(m, i_)
                         * Ginv[bl2](iW_ - iw_)(o, k_);
      }

    return Fpp_iw;
  }

  chi4_iw_t Fph_from_G2ph_conn(chi4_iw_t::const_view_type G2ph_conn_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();

    // Temporary quantities
    g_iw_t Ginv = inverse(G_iw);

    // Calculate vertex function F
    chi4_iw_t Fph_iw = G2ph_conn_iw; // FIXME Product Ranges with += Lazy Expressions
    Fph_iw()         = 0;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = G2ph_conn_iw(bl1, bl2).target_shape()[0];
        int bl2_size = G2ph_conn_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                Fph_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << Fph_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + Ginv[bl1](iW_ + iw_)(j_, n) * Ginv[bl2](iwp_)(l_, p) * G2ph_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p) * Ginv[bl1](iw_)(m, i_)
                         * Ginv[bl2](iW_ + iwp_)(o, k_);
      }

    return Fph_iw;
  }

  chi4_iw_t G2_from_G2c(chi4_iw_t::const_view_type G2_conn_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Calculate G2_iw from G2_conn_iw and G_iw
    chi4_iw_t G2_iw = G2_conn_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        G2_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << G2_conn_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
              + beta * kronecker(iw1_, iw2_) * G_iw[bl1](iw1_)(j_, i_) * G_iw[bl2](iw3_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iw2_, iw3_) * G_iw[bl1](iw1_)(l_, i_) * G_iw[bl2](iw3_)(j_, k_);

    return G2_iw;
  }

  chi4_iw_t G2pp_from_G2pp_conn(chi4_iw_t::const_view_type G2pp_conn_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Calculate G2_iw from G2_conn_iw and G_iw
    chi4_iw_t G2pp_iw = G2pp_conn_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        G2pp_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2pp_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
              + beta * kronecker(iw_, iW_ - iwp_) * G_iw[bl1](iw_)(j_, i_) * G_iw[bl2](iW_ - iw_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iW_ - iwp_, iW_ - iw_) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ - iw_)(j_, k_);

    return G2pp_iw;
  }

  chi4_iw_t G2ph_from_G2ph_conn(chi4_iw_t::const_view_type G2ph_conn_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Calculate G2_iw from G2_conn_iw and G_iw
    chi4_iw_t G2ph_iw = G2ph_conn_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        G2ph_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2ph_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
              + beta * kronecker(iw_, iW_ + iw_) * G_iw[bl1](iw_)(j_, i_) * G_iw[bl2](iW_ + iwp_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iW_ + iw_, iW_ + iwp_) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ + iwp_)(j_, k_);

    return G2ph_iw;
  }

  chi4_iw_t chi_tilde_ph_from_G2ph_conn(chi4_iw_t::const_view_type G2ph_conn_iw, g_iw_cv_t G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Calculate chi_tilde_pha from G2ph_conn_iw and G_iw
    chi4_iw_t chi_tilde_ph = G2ph_conn_iw;

    // Calculate generalized susceptibility in the ph channel from G2_conn_iw and G_iw
    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        chi_tilde_ph(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2ph_conn_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
              - beta * kronecker(bl1, bl2) * kronecker(iw_, iwp_) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iwp_ + iW_)(j_, k_);

    return chi_tilde_ph;
  }

} // namespace triqs_ctint
