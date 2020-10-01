#include "post_process.hpp"
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

namespace triqs_ctint {

  chi4_iw_t G2c_from_M4(chi4_iw_t::const_view_type M4_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Calculate connected part of M4
    chi4_iw_t M4_iw_conn = M4_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        M4_iw_conn(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << M4_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
              - beta * kronecker(iw1_, iw2_) * M_iw[bl1](iw1_)(j_, i_) * M_iw[bl2](iw3_)(l_, k_)
              + beta * kronecker(bl1, bl2) * kronecker(iw2_, iw3_) * M_iw[bl1](iw1_)(l_, i_) * M_iw[bl2](iw3_)(j_, k_);

    // Calculate disconnected part of the two-particle Green function
    chi4_iw_t G2c_iw = M4_iw_conn; // FIXME Product Ranges with += Lazy Expressions
    G2c_iw()         = 0.;

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

  chi4_iw_t G2_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_t::const_view_type G_iw) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].domain().beta;

    // Calculate G2_iw from G2c_iw and G_iw
    chi4_iw_t G2_iw = G2c_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        G2_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << G2c_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
              + beta * kronecker(iw1_, iw2_) * G_iw[bl1](iw1_)(j_, i_) * G_iw[bl2](iw3_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iw2_, iw3_) * G_iw[bl1](iw1_)(l_, i_) * G_iw[bl2](iw3_)(j_, k_);

    return G2_iw;
  }

  chi4_iw_t chi_tilde_ph_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_cv_t G_iw, gf_struct_t const &gf_struct) {

    int n_blocks = G_iw.size();
    double beta  = G_iw[0].domain().beta;

    // Choose ranges such that
    int n_iw_G2      = G2c_iw(0, 0).data().shape()[0] / 2;
    int n_iW         = n_iw_G2 / 2;
    int n_iw         = n_iw_G2 / 2;
    auto imfreq_bos  = gf_mesh<imfreq>{beta, Boson, n_iW};
    auto imfreq_ferm = gf_mesh<imfreq>{beta, Fermion, n_iw};
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
