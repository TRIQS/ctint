#pragma once
#include "./qmc_config.hpp"
#include "./params.hpp"
#include "./fourier.hpp"

namespace triqs_ctint {

  /// Calculate the connected part of the two-particle Green function from M4_iw and M_iw
  chi4_iw_t G2c_from_M4(chi4_iw_t::const_view_type M4_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw);

  /// Calculate the vertex function $F$ from G2c_iw and G_iw
  chi4_iw_t F_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_t::const_view_type G_iw);

  /// Calculate the two-particle Green function from G2c_iw and G_iw
  chi4_iw_t G2_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_t::const_view_type G_iw);

  /// Calculate the $\chi_3$ function from the building blocks M3_iw and M_iw
  template <Chan_t Chan> chi3_iw_t chi3_from_M3(chi3_iw_t::const_view_type M3_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Connected part of M3
    chi3_iw_t M3_iw_conn = M3_iw;

    // Connected part of chi3
    chi3_iw_t chi3_iw = M3_iw;
    chi3_iw()         = 0.;

    // Temporary quantities
    g_iw_t GM   = G0_iw * M_iw;
    g_iw_t MG   = M_iw * G0_iw;
    g_iw_t GMG  = G0_iw * M_iw * G0_iw;
    g_iw_t G_iw = G0_iw + G0_iw * M_iw * G0_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        // Capture block-sizes
        int bl1_size = M3_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M3_iw(bl1, bl2).target_shape()[2];

        if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

          M3_iw_conn(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                - GM[bl1](iw_)(j_, i_) * GM[bl2](iW_ - iw_)(l_, k_) + kronecker(bl1, bl2) * GM[bl1](iw_)(l_, i_) * GM[bl2](iW_ - iw_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl2_size))
              chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                    + G0_iw[bl1](iw_)(m, i_) * G0_iw[bl2](iW_ - iw_)(n, k_) * M3_iw_conn(bl1, bl2)(iw_, iW_)(m, j_, n, l_);

          // Disconnected part
          chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                + G_iw[bl1](iw_)(j_, i_) * G_iw[bl2](iW_ - iw_)(l_, k_) - kronecker(bl1, bl2) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ - iw_)(j_, k_);

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

          M3_iw_conn(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                - beta * kronecker(iW_) * M_iw[bl1](iw_)(j_, i_) * density(GMG[bl2])(l_, k_)
                + kronecker(bl1, bl2) * GM[bl1](iw_)(l_, i_) * MG[bl2](iW_ + iw_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl1_size))
              chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                    + G0_iw[bl1](iw_)(m, i_) * G0_iw[bl1](iW_ + iw_)(j_, n) * M3_iw_conn(bl1, bl2)(iw_, iW_)(m, n, k_, l_);

          // Disconnected part
          chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                + beta * kronecker(iW_) * G_iw[bl1](iw_)(j_, i_) * density(G_iw[bl2])(l_, k_)
                - kronecker(bl1, bl2) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ + iw_)(j_, k_);
        }
      }

    return chi3_iw;
  }

} // namespace triqs_ctint
