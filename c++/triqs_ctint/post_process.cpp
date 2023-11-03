#include "post_process.hpp"
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

namespace triqs_ctint {

  chi4_iw_t G2c_from_M4(chi4_iw_t::const_view_type M4_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    chi4_iw_t G2c_iw = M4_iw; // FIXME Product Ranges with += Lazy Expressions

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
    G2c_iw() = 0.;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = M4_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M4_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                G2c_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << G2c_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
                      + G0_iw[bl1](iw2_)(j_, n) * G0_iw[bl2](iw1_ - iw2_ + iw3_)(l_, p) * M4_iw_conn(bl1, bl2)(iw1_, iw2_, iw3_)(m, n, o, p)
                         * G0_iw[bl1](iw1_)(m, i_) * G0_iw[bl2](iw3_)(o, k_);
      }

    return G2c_iw;
  }

  chi4_iw_t G2ppc_from_M4pp(chi4_iw_t::const_view_type M4pp_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    chi4_iw_t G2ppc_iw = M4pp_iw; // FIXME Product Ranges with += Lazy Expressions

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
    G2ppc_iw() = 0.;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = M4pp_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M4pp_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                G2ppc_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2ppc_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + G0_iw[bl1](iW_ - iwp_)(j_, n) * G0_iw[bl2](iwp_)(l_, p) * M4pp_iw_conn(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p)
                         * G0_iw[bl1](iw_)(m, i_) * G0_iw[bl2](iW_ - iw_)(o, k_);
      }

    return G2ppc_iw;
  }

  chi4_iw_t G2phc_from_M4ph(chi4_iw_t::const_view_type M4ph_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    chi4_iw_t G2phc_iw = M4ph_iw; // FIXME Product Ranges with += Lazy Expressions

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
    G2phc_iw() = 0.;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = M4ph_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M4ph_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                G2phc_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2phc_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + G0_iw[bl1](iW_ + iw_)(j_, n) * G0_iw[bl2](iwp_)(l_, p) * M4ph_iw_conn(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p)
                         * G0_iw[bl1](iw_)(m, i_) * G0_iw[bl2](iW_ + iwp_)(o, k_);
      }

    return G2phc_iw;
  }

  chi4_iw_t F_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();

    // Temporary quantities
    g_iw_t Ginv = inverse(G_iw);

    // Calculate vertex function F
    chi4_iw_t F_iw = G2c_iw; // FIXME Product Ranges with += Lazy Expressions
    F_iw()         = 0;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = G2c_iw(bl1, bl2).target_shape()[0];
        int bl2_size = G2c_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                F_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << F_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
                      + Ginv[bl1](iw2_)(j_, n) * Ginv[bl2](iw1_ - iw2_ + iw3_)(l_, p) * G2c_iw(bl1, bl2)(iw1_, iw2_, iw3_)(m, n, o, p)
                         * Ginv[bl1](iw1_)(m, i_) * Ginv[bl2](iw3_)(o, k_);
      }

    return F_iw;
  }

  chi4_iw_t Fpp_from_G2ppc(chi4_iw_t::const_view_type G2ppc_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();

    // Temporary quantities
    g_iw_t Ginv = inverse(G_iw);

    // Calculate vertex function F
    chi4_iw_t Fpp_iw = G2ppc_iw; // FIXME Product Ranges with += Lazy Expressions
    Fpp_iw()         = 0;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = G2ppc_iw(bl1, bl2).target_shape()[0];
        int bl2_size = G2ppc_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                Fpp_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << Fpp_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + Ginv[bl1](iW_ - iwp_)(j_, n) * Ginv[bl2](iwp_)(l_, p) * G2ppc_iw(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p) * Ginv[bl1](iw_)(m, i_)
                         * Ginv[bl2](iW_ - iw_)(o, k_);
      }

    return Fpp_iw;
  }

  chi4_iw_t Fph_from_G2phc(chi4_iw_t::const_view_type G2phc_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();

    // Temporary quantities
    g_iw_t Ginv = inverse(G_iw);

    // Calculate vertex function F
    chi4_iw_t Fph_iw = G2phc_iw; // FIXME Product Ranges with += Lazy Expressions
    Fph_iw()         = 0;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = G2phc_iw(bl1, bl2).target_shape()[0];
        int bl2_size = G2phc_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                Fph_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << Fph_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
                      + Ginv[bl1](iW_ + iw_)(j_, n) * Ginv[bl2](iwp_)(l_, p) * G2phc_iw(bl1, bl2)(iW_, iw_, iwp_)(m, n, o, p) * Ginv[bl1](iw_)(m, i_)
                         * Ginv[bl2](iW_ + iwp_)(o, k_);
      }

    return Fph_iw;
  }

  chi4_iw_t Fpp_loc_from_f4pp_loc(chi4_iw_t::const_view_type f4pp_loc_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type GinvG01_iw,
                                  g_iw_t::const_view_type GinvG02_iw) {

    // calculate Ginv * G0 * M * G0 * Ginv
    auto m = g_iw_t{M_iw};
    m()    = 0;

    for (int bl : range(M_iw.size())) {
      int bl_size      = m[bl].target_shape()[0];
      auto &m_bl       = m[bl];
      auto const &M_bl = M_iw[bl];
      auto const &L_bl = GinvG01_iw[bl];
      auto const &R_bl = GinvG02_iw[bl];

      for (auto iw : m_bl.mesh())
        for (int i : range(bl_size))
          for (int j : range(bl_size))
            for (int k : range(bl_size)) { m_bl[iw](i, i) += L_bl(iw)(i, j) * M_bl[iw](j, k) * R_bl(iw)(k, i); }
    }

    // calculate Fpp_loc
    auto Fpp_loc_iw     = chi4_iw_t{f4pp_loc_iw};
    auto const &iW_mesh = std::get<0>(Fpp_loc_iw(0, 0).mesh());
    auto const &iw_mesh = std::get<1>(Fpp_loc_iw(0, 0).mesh());
    auto const beta     = iW_mesh.beta();

    for (int bl1 : range(M_iw.size()))
      for (int bl2 : range(M_iw.size())) {
        int bl_size    = m[bl1].target_shape()[0];
        auto const &m1 = m[bl1];
        auto const &m2 = m[bl2];
        auto &F_loc    = Fpp_loc_iw(bl1, bl2);

        for (auto iW : iW_mesh)
          for (auto iw : iw_mesh)
            for (auto iwp : iw_mesh)
              for (int i : range(bl_size)) {
                F_loc[iW, iw, iwp](i, i, i, i) -= beta * kronecker(iw, iW - iwp) * m1(iw)(i, i) * m2(iW - iw)(i, i);
                if (bl1 == bl2) { F_loc[iW, iw, iwp](i, i, i, i) += beta * kronecker(iW - iwp, iW - iw) * m1(iw)(i, i) * m2(iW - iw)(i, i); }
              }
      }

    return Fpp_loc_iw;
  }

  chi4_iw_t Fph_loc_from_f4ph_loc(chi4_iw_t::const_view_type f4ph_loc_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type GinvG01_iw,
                                  g_iw_t::const_view_type GinvG02_iw) {

    // calculate Ginv * G0 * M * G0 * Ginv
    auto m = g_iw_t{M_iw};
    m()    = 0;

    for (int bl : range(M_iw.size())) {
      int bl_size      = m[bl].target_shape()[0];
      auto &m_bl       = m[bl];
      auto const &M_bl = M_iw[bl];
      auto const &L_bl = GinvG01_iw[bl];
      auto const &R_bl = GinvG02_iw[bl];

      for (auto iw : m_bl.mesh())
        for (int i : range(bl_size))
          for (int j : range(bl_size))
            for (int k : range(bl_size)) { m_bl[iw](i, i) += L_bl(iw)(i, j) * M_bl[iw](j, k) * R_bl(iw)(k, i); }
    }

    // calculate Fph_loc
    auto Fph_loc_iw     = chi4_iw_t{f4ph_loc_iw};
    auto const &iW_mesh = std::get<0>(Fph_loc_iw(0, 0).mesh());
    auto const &iw_mesh = std::get<1>(Fph_loc_iw(0, 0).mesh());
    auto const beta     = iW_mesh.beta();

    for (int bl1 : range(M_iw.size()))
      for (int bl2 : range(M_iw.size())) {
        int bl_size    = m[bl1].target_shape()[0];
        auto const &m1 = m[bl1];
        auto const &m2 = m[bl2];
        auto &F_loc    = Fph_loc_iw(bl1, bl2);

        for (auto iW : iW_mesh)
          for (auto iw : iw_mesh)
            for (auto iwp : iw_mesh)
              for (int i : range(bl_size)) {
                F_loc[iW, iw, iwp](i, i, i, i) -= beta * kronecker(iw, iW + iw) * m1(iw)(i, i) * m2(iW + iwp)(i, i);
                if (bl1 == bl2) { F_loc[iW, iw, iwp](i, i, i, i) += beta * kronecker(iW + iw, iW + iwp) * m1(iw)(i, i) * m2(iW + iwp)(i, i); }
              }
      }

    return Fph_loc_iw;
  }

  chi4_iw_t G2_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Calculate G2_iw from G2c_iw and G_iw
    chi4_iw_t G2_iw = G2c_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        G2_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << G2c_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
              + beta * kronecker(iw1_, iw2_) * G_iw[bl1](iw1_)(j_, i_) * G_iw[bl2](iw3_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iw2_, iw3_) * G_iw[bl1](iw1_)(l_, i_) * G_iw[bl2](iw3_)(j_, k_);

    return G2_iw;
  }

  chi4_iw_t G2pp_from_G2ppc(chi4_iw_t::const_view_type G2ppc_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Calculate G2_iw from G2c_iw and G_iw
    chi4_iw_t G2pp_iw = G2ppc_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        G2pp_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2ppc_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
              + beta * kronecker(iw_, iW_ - iwp_) * G_iw[bl1](iw_)(j_, i_) * G_iw[bl2](iW_ - iw_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iW_ - iwp_, iW_ - iw_) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ - iw_)(j_, k_);

    return G2pp_iw;
  }

  chi4_iw_t G2ph_from_G2phc(chi4_iw_t::const_view_type G2phc_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Calculate G2_iw from G2c_iw and G_iw
    chi4_iw_t G2ph_iw = G2phc_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        G2ph_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2phc_iw(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_)
              + beta * kronecker(iw_, iW_ + iw_) * G_iw[bl1](iw_)(j_, i_) * G_iw[bl2](iW_ + iwp_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iW_ + iw_, iW_ + iwp_) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ + iwp_)(j_, k_);

    return G2ph_iw;
  }

  chi4_iw_t chi_tilde_ph_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_cv_t G_iw, gf_struct_t const &gf_struct, int n_iW = -1, int n_iw = -1) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].mesh().beta();

    // Choose ranges such that
    int n_iw_G2 = G2c_iw(0, 0).data().shape()[0] / 2;
    n_iW        = (n_iW < 0 ? n_iw_G2 / 2 : n_iW);
    n_iw        = (n_iw < 0 ? n_iw_G2 / 2 : n_iw);
    TRIQS_ASSERT(n_iw + n_iW <= n_iw_G2);
    auto imfreq_bos  = mesh::imfreq{beta, Boson, n_iW};
    auto imfreq_ferm = mesh::imfreq{beta, Fermion, n_iw};
    auto mesh        = prod{imfreq_bos, imfreq_ferm, imfreq_ferm};

    chi4_iw_t chi_tilde_ph = make_block2_gf(mesh, gf_struct);

    // Calculate generalized susceptibility in the ph channel from G2c_iw and G_iw
    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        chi_tilde_ph(bl1, bl2)(iW_, iw_, iwp_)(i_, j_, k_, l_) << G2c_iw(bl1, bl2)(iw_, iw_ + iW_, iwp_ + iW_)(i_, j_, k_, l_)
              - beta * kronecker(bl1, bl2) * kronecker(iw_, iwp_) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iwp_ + iW_)(j_, k_);

    return chi_tilde_ph;
  }

} // namespace triqs_ctint
