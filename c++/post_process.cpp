#include "post_process.hpp"
#include <triqs/gfs.hpp>
#include <triqs/gfs/singularity/fit_tail.hpp>

//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

namespace triqs_ctint {

  chi4_iw_t F_from_M4(chi4_iw_t::const_view_type M4_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    g_iw_t G_iw = G0_iw + G0_iw * M_iw * G0_iw;
    double beta = M_iw[0].domain().beta;

    // The connected part of M4
    chi4_iw_t M4_iw_conn = M4_iw;
    M4_iw_conn(bl1_, bl2_)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << M4_iw(bl1_, bl2_)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
          + beta * kronecker(iw1_, iw2_) * M_iw[bl1_](iw1_)(j_, i_) * M_iw[bl2_](iw3_)(l_, k_)
          - beta * kronecker(bl1_, bl2_) * kronecker(iw2_, iw3_) * M_iw[bl1_](iw1_)(l_, i_) * M_iw[bl2_](iw3_)(j_, k_);

    // Temporary quantities
    auto Ginv_x_G0 = inverse(G_iw) * G0_iw;
    auto G0_x_Ginv = G0_iw * inverse(G_iw);

    // The vertex function F
    chi4_iw_t F_iw = M4_iw_conn; // FIXME Product Ranges with += Lazy Expressions and or Einstein Summation
    //F_iw(bl1_, bl2_)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << Ginv_x_G0[bl1_](iw1_)(j_, n_) * Ginv_x_G0[bl2_](iw1_ + iw2_ - iw3_)(l_, p_)
    //* M4_iw_conn(bl1_, bl2_)(iw1_, iw2_, iw3_)(m_, n_, o_, p_) * G0_x_Ginv[bl1_](iw1_)(m_, i_) * G0_x_Ginv[bl2_](iw3_)(o_, k_);

    return F_iw;
  }

  chi3_iw_t K2_from_M3(chi3_iw_t::const_view_type M3_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    chi3_iw_t K2_iw = M3_iw;
    //K2_iw(bl1_, bl2_)(iw1_, iw3_)(i_, j_, k_, l_) << Ginv_x_G0[bl1_](iw1_)(j_, n_) * Ginv_x_G0[bl2_](iw1_ + iw2_ - iw3_)(l_, p_)
    //* M4_iw_conn(bl1_, bl2_)(iw1_, iw2_, iw3_)(m_, n_, o_, p_) * G0_x_Ginv[bl1_](iw1_)(m_, i_) * G0_x_Ginv[bl2_](iw3_)(o_, k_);

    return K2_iw;
  }

} // namespace triqs_ctint
