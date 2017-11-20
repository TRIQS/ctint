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

  // Calculate the $\chi_2$ function from the building blocks M2_tau and M_iw
  template <Chan_t Chan> chi2_tau_t chi2_from_M2(chi2_tau_t::const_view_type M2_tau, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Chi2 and its connected part
    chi2_tau_t chi2_tau      = M2_tau;
    chi2_tau_t chi2_tau_conn = M2_tau;

    // Temporary quantities
    g_iw_t G_iw = G0_iw + G0_iw * M_iw * G0_iw;

#ifdef GTAU_IS_COMPLEX
    auto GMG_tau  = make_gf_from_inverse_fourier(g_iw_t{G0_iw * M_iw * G0_iw}, 10000);
    g_tau_t G_tau = make_gf_from_inverse_fourier(G_iw, 10000);
#else
    auto GMG_tau  = get_real(make_gf_from_inverse_fourier(g_iw_t{G0_iw * M_iw * G0_iw}, 10000), true);
    g_tau_t G_tau = get_real(make_gf_from_inverse_fourier(G_iw, 10000), true);
#endif

    auto const &tau_mesh = M2_tau(0, 0).mesh();

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        // Capture block-sizes
        int bl1_size = M2_tau(bl1, bl2).target_shape()[0];
        int bl2_size = M2_tau(bl1, bl2).target_shape()[2];

        if (Chan == Chan_t::PP) { // =====  Particle-particle channel // FIXME c++17 if constexpr

          //chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_) << M2_tau(bl1, bl2)(t_)(i_, j_, k_, l_) // FIXME TAIL ERROR
          //- GMG_tau[bl1](beta - t_)(j_, i_) * GMG_tau[bl2](beta - t_)(l_, k_)
          //+ kronecker(bl1, bl2) * GMG_tau[bl1](beta - t_)(l_, i_) * GMG_tau[bl2](beta - t_)(j_, k_);

          //chi2_tau(bl1, bl2)(t_)(i_, j_, k_, l_) << chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_)
          //+ G_tau[bl1](beta - t_)(j_, i_) * G_tau[bl2](beta - t_)(l_, k_)
          //+ kronecker(bl1, bl2) * G_tau[bl1](beta - t_)(l_, i_) * G_tau[bl2](beta - t_)(j_, k_);

          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size))
                  for (auto const &t : tau_mesh)
                    chi2_tau_conn(bl1, bl2)[t](i, j, k, l) = M2_tau(bl1, bl2)[t](i, j, k, l)
                       - GMG_tau[bl1](beta - t)(j, i) * GMG_tau[bl2](beta - t)(l, k)
                       + kronecker(bl1, bl2) * GMG_tau[bl1](beta - t)(l, i) * GMG_tau[bl2](beta - t)(j, k);

          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size))
                  for (auto const &t : tau_mesh)
                    chi2_tau(bl1, bl2)[t](i, j, k, l) = chi2_tau_conn(bl1, bl2)[t](i, j, k, l)
                       + G_tau[bl1](beta - t)(j, i) * G_tau[bl2](beta - t)(l, k) // Disconnected part
                       - kronecker(bl1, bl2) * G_tau[bl1](beta - t)(l, i) * G_tau[bl2](beta - t)(j, k);

        } else if (Chan == Chan_t::PH) { // ===== Particle-hole channel

          //chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_) << M2_tau(bl1, bl2)(t_)(i_, j_, k_, l_)
          //- GMG_tau[bl1](beta - 1e-10)(j_, i_) * GMG_tau[bl2](beta - 1e-10)(l_, k_)
          //- kronecker(bl1, bl2) * GMG_tau[bl1](beta - t_)(l_, i_) * GMG_tau[bl2](t_)(j_, k_);

          //chi2_tau(bl1, bl2)(t_)(i_, j_, k_, l_) << chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_)
          //+ G_tau[bl1](beta - 1e-10)(j_, i_) * G_tau[bl2](beta - 1e-10)(l_, k_)
          //+ kronecker(bl1, bl2) * G_tau[bl1](beta - t_)(l_, i_) * G_tau[bl2](t_)(j_, k_);

          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size))
                  for (auto const &t : tau_mesh)
                    chi2_tau_conn(bl1, bl2)[t](i, j, k, l) = M2_tau(bl1, bl2)[t](i, j, k, l)
                       - GMG_tau[bl1](beta - 1e-10)(j, i) * GMG_tau[bl2](beta - 1e-10)(l, k)
                       - kronecker(bl1, bl2) * GMG_tau[bl1](beta - t)(l, i) * GMG_tau[bl2](t)(j, k);

          for (int i : range(bl1_size))
            for (int j : range(bl1_size))
              for (int k : range(bl2_size))
                for (int l : range(bl2_size))
                  for (auto const &t : tau_mesh)
                    chi2_tau(bl1, bl2)[t](i, j, k, l) = chi2_tau_conn(bl1, bl2)[t](i, j, k, l)
                       + G_tau[bl1](beta - 1e-10)(j, i) * G_tau[bl2](beta - 1e-10)(l, k) // Disconnected part
                       + kronecker(bl1, bl2) * G_tau[bl1](beta - t)(l, i) * G_tau[bl2](t)(j, k);
        }
      }

    return chi2_tau;
  }

  /// Calculate the $\chi_3$ function from the building blocks M3_iw and M_iw
  template <Chan_t Chan> chi3_iw_t chi3_from_M3(chi3_iw_t::const_view_type M3_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw) {

    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Connected part of M3
    chi3_iw_t M3_iw_conn = M3_iw;

    // Connected part of chi3
    chi3_iw_t chi3_iw = M3_iw;

    // Temporary quantities
    auto G0_x_M = G0_iw * M_iw;
    auto M_x_G0 = M_iw * G0_iw;
    g_iw_t G_iw = G0_iw + G0_iw * M_iw * G0_iw;

#ifdef GTAU_IS_COMPLEX
    auto GMG_tau  = make_gf_from_inverse_fourier(g_iw_t{G0_iw * M_iw * G0_iw}, 10000);
    g_tau_t G_tau = make_gf_from_inverse_fourier(G_iw, 10000);
#else
    auto GMG_tau  = get_real(make_gf_from_inverse_fourier(g_iw_t{G0_iw * M_iw * G0_iw}, 10000), true);
    g_tau_t G_tau = get_real(make_gf_from_inverse_fourier(G_iw, 10000), true);
#endif

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        // Capture block-sizes
        int bl1_size = M3_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M3_iw(bl1, bl2).target_shape()[2];

        if (Chan == Chan_t::PP) { // =====  Particle-particle channel // FIXME c++17 if constexpr

          M3_iw_conn(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                - G0_x_M[bl1](iw_)(j_, i_) * G0_x_M[bl2](iW_ - iw_)(l_, k_) + kronecker(bl1, bl2) * G0_x_M[bl1](iw_)(l_, i_) * G0_x_M[bl2](iW_ - iw_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl2_size))
              chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                 << G0_iw[bl1](iw_)(m, i_) * G0_iw[bl2](iW_ - iw_)(n, k_) * M3_iw_conn(bl1, bl2)(iw_, iW_)(m, j_, n, l_)
                    + G_iw[bl1](iw_)(j_, i_) * G_iw[bl2](iW_ - iw_)(l_, k_) // Disconnected part
                    - kronecker(bl1, bl2) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ - iw_)(j_, k_);

        } else if (Chan == Chan_t::PH) { // ===== Particle-hole channel

          M3_iw_conn(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                + beta * kronecker(iw_, iW_ + iw_) * M_iw[bl1](iw_)(j_, i_) * GMG_tau[bl2](beta - 1e-10)(l_, k_)
                + kronecker(bl1, bl2) * G0_x_M[bl1](iw_)(l_, i_) * M_x_G0[bl2](iW_ + iw_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl1_size))
              chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_)
                 << G0_iw[bl1](iw_)(m, i_) * G0_iw[bl1](iW_ + iw_)(j_, n) * M3_iw_conn(bl1, bl2)(iw_, iW_)(m, n, k_, l_)
                    - beta * kronecker(iw_, iW_ + iw_) * G_iw[bl1](iw_)(j_, i_) * G_tau[bl2](beta - 1e-10)(l_, k_) // Disconnected part
                    - kronecker(bl1, bl2) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ + iw_)(j_, k_);
        }
      }

    return chi3_iw;
  }

} // namespace triqs_ctint
