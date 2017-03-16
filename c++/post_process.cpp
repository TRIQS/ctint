#include "post_process.hpp"
#include <triqs/gfs.hpp>
#include <triqs/gfs/singularity/fit_tail.hpp>

//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

namespace triqs_ctint {

  chi4_iw_t F_from_M4(chi4_iw_t::const_view_type M4_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    g_iw_t G_iw  = G0_iw + G0_iw * M_iw * G0_iw;
    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Calculate connected part of M4
    chi4_iw_t M4_iw_conn = M4_iw;
    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks))
        M4_iw_conn(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << M4_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
              + beta * kronecker(iw1_, iw2_) * M_iw[bl1](iw1_)(j_, i_) * M_iw[bl2](iw3_)(l_, k_)
              - beta * kronecker(bl1, bl2) * kronecker(iw2_, iw3_) * M_iw[bl1](iw1_)(l_, i_) * M_iw[bl2](iw3_)(j_, k_);

    // Temporary quantities
    auto Ginv_x_G0 = inverse(G_iw) * G0_iw;
    auto G0_x_Ginv = G0_iw * inverse(G_iw);

    // Calculate vertex function F
    chi4_iw_t F_iw = M4_iw_conn; // FIXME Product Ranges with += Lazy Expressions
    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        int bl1_size = M4_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M4_iw(bl1, bl2).target_shape()[2];

        for (int m : range(bl1_size))
          for (int n : range(bl1_size))
            for (int o : range(bl2_size))
              for (int p : range(bl2_size))
                F_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << F_iw(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
                      + Ginv_x_G0[bl1](iw1_)(j_, n) * Ginv_x_G0[bl2](iw1_ + iw2_ - iw3_)(l_, p) * M4_iw_conn(bl1, bl2)(iw1_, iw2_, iw3_)(m, n, o, p)
                         * G0_x_Ginv[bl1](iw1_)(m, i_) * G0_x_Ginv[bl2](iw3_)(o, k_);
      }

    return F_iw;
  }

} // namespace triqs_ctint
