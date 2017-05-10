#pragma once
#include "./qmc_config.hpp"
#include "./params.hpp"
#include "./fourier.hpp"

namespace triqs_ctint {

  /// Calculate the vertex function $F$ from the the building blocks M4_iw and M_iw
  chi4_iw_t F_from_M4(chi4_iw_t::const_view_type M4_iw, block_gf_const_view<imfreq, matrix_valued> M_iw,
                      block_gf_const_view<imfreq, matrix_valued> G0_iw);

  /// Calculate the connected part of the two-particle Green function from M4_iw and M_iw
  chi4_iw_t G2c_from_M4(chi4_iw_t::const_view_type M4_iw, block_gf_const_view<imfreq, matrix_valued> M_iw,
                      block_gf_const_view<imfreq, matrix_valued> G0_iw);


  /// Calculate the $\chi_3$ function from the building block M3_iw and M_iw
  template <Chan_t Chan> chi3_iw_t chi3_from_M3(chi3_iw_t::const_view_type M3_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Connected part of M3
    chi3_iw_t M3_iw_conn = M3_iw;

    // Connected part of chi3
    chi3_iw_t chi3_iw = M3_iw;

    // Temporary quantities
    auto G0_x_M  = G0_iw * M_iw;
    auto M_x_G0  = M_iw * G0_iw;
    auto GMG_tau = make_gf_from_inverse_fourier(g_iw_t{G0_iw * M_iw * G0_iw}, 10000);

    g_iw_t G_iw   = G0_iw + G0_iw * M_iw * G0_iw;
    g_tau_t G_tau = make_gf_from_inverse_fourier(G_iw, 10000);

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        // Capture block-sizes
        int bl1_size = M3_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M3_iw(bl1, bl2).target_shape()[2];

        if (Chan == Chan_t::PP) { // =====  Particle-particle channel // FIXME c++17 if constexpr

          M3_iw_conn(bl1, bl2)(iw1_, iw3_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)(iw1_, iw3_)(i_, j_, k_, l_)
                - G0_x_M[bl1](iw1_)(j_, i_) * G0_x_M[bl2](iw3_)(l_, k_) + kronecker(bl1, bl2) * G0_x_M[bl1](iw1_)(l_, i_) * G0_x_M[bl2](iw3_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl2_size))
              chi3_iw(bl1, bl2)(iw1_, iw3_)(i_, j_, k_, l_)
                 << G0_iw[bl1](iw1_)(m, i_) * G0_iw[bl2](iw3_)(n, k_) * M3_iw_conn(bl1, bl2)(iw1_, iw3_)(m, j_, n, l_)
                    + G_iw[bl1](iw1_)(j_, i_) * G_iw[bl2](iw3_)(l_, k_)
                    - kronecker(bl1, bl2) * G_iw[bl1](iw1_)(l_, i_) * G_iw[bl2](iw3_)(j_, k_); // Disconnected part

        } else if (Chan == Chan_t::PH) { // ===== Particle-hole channel

          M3_iw_conn(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) // FIXME Check 0-
                - beta * kronecker(iw1_, iw2_) * M_iw[bl1](iw1_)(j_, i_) * GMG_tau[bl2](-1e-10)(l_, k_)
                + kronecker(bl1, bl2) * G0_x_M[bl1](iw1_)(l_, i_) * M_x_G0[bl2](iw2_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl1_size))
              chi3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_)
                 << G0_iw[bl1](iw1_)(m, i_) * G0_iw[bl1](iw2_)(j_, n) * M3_iw_conn(bl1, bl2)(iw1_, iw2_)(m, n, k_, l_)
                    + beta * kronecker(iw1_, iw2_) * G_iw[bl1](iw1_)(j_, i_) * G_tau[bl2](-1e-10)(l_, k_) // Disconnected part
                    - kronecker(bl1, bl2) * G_iw[bl1](iw1_)(l_, i_) * G_iw[bl2](iw2_)(j_, k_);
        }
      }

    return chi3_iw;
  }

} // namespace triqs_ctint
