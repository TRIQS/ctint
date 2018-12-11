#pragma once
#include "./qmc_config.hpp"
#include "./params.hpp"

namespace triqs_ctint {

  /// Calculate the connected part of the two-particle Green function from M4_iw and M_iw
  chi4_iw_t G2c_from_M4(chi4_iw_t::const_view_type M4_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw);

  /// Calculate the vertex function $F$ from G2c_iw and G_iw
  chi4_iw_t F_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_cv_t G_iw);

  /// Calculate the two-particle Green function from G2c_iw and G_iw
  chi4_iw_t G2_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_cv_t G_iw);

  /// Calculate the $\chi_3$ function from the building blocks M3_iw and M_iw
  template <Chan_t Chan> chi3_iw_t chi3_from_M3(chi3_iw_t::const_view_type M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) {

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

          M3_iw_conn(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)[iw_, iW_](i_, j_, k_, l_)
                - GM[bl1](iw_)(j_, i_) * GM[bl2](iW_ - iw_)(l_, k_) + kronecker(bl1, bl2) * GM[bl1](iw_)(l_, i_) * GM[bl2](iW_ - iw_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl2_size))
              chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw_, iW_](i_, j_, k_, l_)
                    + G0_iw[bl1](iw_)(m, i_) * G0_iw[bl2](iW_ - iw_)(n, k_) * M3_iw_conn(bl1, bl2)(iw_, iW_)(m, j_, n, l_);

          // Disconnected part
          chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw_, iW_](i_, j_, k_, l_)
                + G_iw[bl1](iw_)(j_, i_) * G_iw[bl2](iW_ - iw_)(l_, k_) - kronecker(bl1, bl2) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ - iw_)(j_, k_);

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

          auto km = make_zero_tail(GMG[bl2], 2);
          M3_iw_conn(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)[iw_, iW_](i_, j_, k_, l_)
                - beta * kronecker(iW_) * M_iw[bl1](iw_)(j_, i_) * density(GMG[bl2], km)(l_, k_)
                + kronecker(bl1, bl2) * GM[bl1](iw_)(l_, i_) * MG[bl2](iW_ + iw_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl1_size))
              chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw_, iW_](i_, j_, k_, l_)
                    + G0_iw[bl1](iw_)(m, i_) * G0_iw[bl1](iW_ + iw_)(j_, n) * M3_iw_conn(bl1, bl2)(iw_, iW_)(m, n, k_, l_);

          // Disconnected part
          matrix_view<dcomplex>{km(1, ellipsis())} = 1.0;
          chi3_iw(bl1, bl2)(iw_, iW_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw_, iW_](i_, j_, k_, l_)
                + beta * kronecker(iW_) * G_iw[bl1](iw_)(j_, i_) * density(G_iw[bl2], km)(l_, k_)
                - kronecker(bl1, bl2) * G_iw[bl1](iw_)(l_, i_) * G_iw[bl2](iW_ + iw_)(j_, k_);
        }
      }

    return chi3_iw;
  }

  // Calculate the $\chi_2$ function from the building blocks M2_tau and M_iw
  template <Chan_t Chan> chi2_tau_t chi2_from_M2(chi2_tau_t::const_view_type M2_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) {

    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Chi2 and its connected part
    chi2_tau_t chi2_tau      = M2_tau;
    chi2_tau_t chi2_tau_conn = M2_tau;

    // Temporary quantities
    g_iw_t GMG_iw = G0_iw * M_iw * G0_iw;
    g_iw_t G_iw   = G0_iw + GMG_iw;

    auto tau_mesh = make_adjoint_mesh(G_iw(0).mesh());

    // TODO
    //auto tail_GMG = make_zero_tail(GMG_iw, 3);
    //for(auto [bl, M_hf_bl]) bl[2] = M_hf_bl;

#ifdef GTAU_IS_COMPLEX
    auto GMG_tau  = make_gf_from_fourier(GMG_iw, tau_mesh, make_zero_tail(GMG_iw, 2));
    g_tau_t G_tau = make_gf_from_fourier(G_iw);
#else
    auto GMG_tau  = real(make_gf_from_fourier(GMG_iw, tau_mesh, make_zero_tail(GMG_iw, 2)));
    g_tau_t G_tau = real(make_gf_from_fourier(G_iw));
    if (!is_gf_real_in_tau(G_iw, 1e-8)) std::cerr << "WARNING: Assuming real G_tau, but found Imag(G(tau)) > 1e-8. Casting to Real.\n";
#endif

    auto dens_GMG = density(GMG_iw, make_zero_tail(GMG_iw, 2));

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

	  //chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_) << M2_tau(bl1, bl2)[t_](i_, j_, k_, l_)
	  //- GMG_tau[bl1](beta - t_)(j_, i_) * GMG_tau[bl2](beta - t_)(l_, k_)
	  //+ kronecker(bl1, bl2) * GMG_tau[bl1](beta - t_)(l_, i_) * GMG_tau[bl2](beta - t_)(j_, k_);

          // Add Disconnected part
          chi2_tau(bl1, bl2)(t_)(i_, j_, k_, l_) << chi2_tau_conn(bl1, bl2)[t_](i_, j_, k_, l_)
                + G_tau[bl1](beta - t_)(j_, i_) * G_tau[bl2](beta - t_)(l_, k_)
                - kronecker(bl1, bl2) * G_tau[bl1](beta - t_)(l_, i_) * G_tau[bl2](beta - t_)(j_, k_);

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

	  //chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_) << M2_tau(bl1, bl2)[t_](i_, j_, k_, l_) - dens_GMG[bl1](j_, i_) * dens_GMG[bl2](l_, k_)
	  //- kronecker(bl1, bl2) * GMG_tau[bl1](beta - t_)(l_, i_) * GMG_tau[bl2](t_)(j_, k_);

          // Add Disconnected part, FIXME used measured density?
          chi2_tau(bl1, bl2)(t_)(i_, j_, k_, l_) << chi2_tau_conn(bl1, bl2)[t_](i_, j_, k_, l_)
                + density(G_iw[bl1])(j_, i_) * density(G_iw[bl2])(l_, k_)
                + kronecker(bl1, bl2) * G_tau[bl1](beta - t_)(l_, i_) * G_tau[bl2](t_)(j_, k_);
        }
      }

    return chi2_tau;
  }

  /// Calculate M3_conn from M3
  template <Chan_t Chan>
  chi3_tau_t M3_conn_from_M3(chi3_tau_t M3_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau, std::vector<matrix<M_tau_scalar_t>> M_hartree) {

    double beta  = M_tau[0].domain().beta;
    int n_blocks = M_tau.size();

    auto const &tau_mesh_M3 = M3_tau(0, 0).mesh();
    double dtau_M3          = std::get<0>(tau_mesh_M3).delta();

    // Temporary quantities
    g_iw_t GM_iw  = G0_iw * M_iw;
    g_iw_t MG_iw  = M_iw * G0_iw;
    g_iw_t GMG_iw = G0_iw * M_iw * G0_iw;

    auto GM = make_gf_from_fourier(GM_iw);
    auto MG = make_gf_from_fourier(MG_iw);

    // Connected part of M3
    chi3_tau_t M3_tau_conn = M3_tau;
    auto dens_GMG          = density(GMG_iw, make_zero_tail(GMG_iw, 2));
    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

          M3_tau_conn(bl1, bl2)(t1_, t2_)(i_, j_, k_, l_) << M3_tau(bl1, bl2)[t1_, t2_](i_, j_, k_, l_)
                - GM[bl1](beta - t1_)(j_, i_) * GM[bl2](beta - t2_)(l_, k_)
                + kronecker(bl1, bl2) * GM[bl1](beta - t1_)(l_, i_) * GM[bl2](beta - t2_)(j_, k_);

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

          for (auto [t1, t2] : tau_mesh_M3) {
            auto [s, d_t2_t1] = cyclic_difference(t2, t1);

            M3_tau_conn(bl1, bl2)[t1, t2](i_, j_, k_, l_) << M3_tau(bl1, bl2)[t1, t2](i_, j_, k_, l_)
                  - (s * M_tau[bl1](d_t2_t1)(j_, i_) + kronecker(t1.index(), t2.index()) / dtau_M3 * M_hartree[bl1](j_, i_)) * dens_GMG[bl2](l_, k_)
                  - kronecker(bl1, bl2) * GM[bl1](beta - t1)(l_, i_) * MG[bl2](t2)(j_, k_); // Sign change from GM shift
          }
        }
      }

    return M3_tau_conn;
  }

  /// Calculate the chi2_tau from M3_tau and M_iw
  template <Chan_t Chan>
  chi2_tau_t M2_from_M3(chi3_tau_t M3_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau, std::vector<matrix<M_tau_scalar_t>> M_hartree,
                        g_tau_cv_t G0_tau, int n_tau_M2) {

    double beta  = G0_tau[0].domain().beta;
    int n_blocks = G0_tau.size();

    auto const &tau_mesh_M3 = M3_tau(0, 0).mesh();
    double dtau_M3          = std::get<0>(tau_mesh_M3).delta();
    int n_tau_M3            = std::get<0>(tau_mesh_M3).size();

    auto tau_mesh_M2 = gf_mesh<imtime>{beta, Boson, n_tau_M2};
    auto M2_tau      = make_block2_gf(tau_mesh_M2, make_const_view(M3_tau));
    M2_tau()         = 0.0;

    chi3_tau_t M3_tau_conn = M3_conn_from_M3<Chan>(M3_tau, M_iw, G0_iw, M_tau, M_hartree);
    //chi3_tau_t M3_tau_conn = M3_tau;

    ////Account for M3 edge bins beeing smaller // FIXME Why not necessary?
    //auto _ = all_t{};
    //int n  = n_tau_M3 - 1;
    //for (auto &M : M3_tau_conn) {
    //M[0, _] *= 0.5;
    //M[_, 0] *= 0.5;
    //M[n, _] *= 0.5;
    //M[_, n] *= 0.5;
    //}

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        // Capture block-sizes
        int bl1_size = M3_tau_conn(bl1, bl2).target_shape()[0];
        int bl2_size = M3_tau_conn(bl1, bl2).target_shape()[2];

        for (auto &t : tau_mesh_M2) {
          auto M2 = M2_tau(bl1, bl2)[t];

          for (auto [t1, t2] : tau_mesh_M3) {
            auto M3 = M3_tau_conn(bl1, bl2)[t1, t2];

            auto [s1, d_t1_t] = cyclic_difference(t1, t);
            auto [s2, d_t2_t] = cyclic_difference(t2, t);
            auto [s3, d_t_t2] = cyclic_difference(t, t2);

            if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

              for (int m : range(bl1_size))
                for (int n : range(bl2_size))
                  M2(i_, j_, k_, l_) << M2(i_, j_, k_, l_) + M3(m, j_, n, l_) * s1 * G0_tau[bl1](d_t1_t)(m, i_) * s2 * G0_tau[bl2](d_t2_t)(n, k_);

            } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel
              for (int m : range(bl1_size))
                for (int n : range(bl1_size))
                  M2(i_, j_, k_, l_) << M2(i_, j_, k_, l_) + M3(m, n, k_, l_) * s1 * G0_tau[bl1](d_t1_t)(m, i_) * s3 * G0_tau[bl1](d_t_t2)(j_, n);
            }
          }
        }
      }

    M2_tau() = M2_tau * dtau_M3 * dtau_M3;
    return M2_tau;
  }

  // For wrapping purposes
  inline chi3_iw_t chi3_from_M3_PP(chi3_iw_t::const_view_type M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) {
    return chi3_from_M3<Chan_t::PP>(M3_iw, M_iw, G0_iw);
  }
  inline chi3_iw_t chi3_from_M3_PH(chi3_iw_t::const_view_type M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) {
    return chi3_from_M3<Chan_t::PH>(M3_iw, M_iw, G0_iw);
  }
  inline chi2_tau_t chi2_from_M2_PP(chi2_tau_t::const_view_type M2_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) {
    return chi2_from_M2<Chan_t::PP>(M2_tau, M_iw, G0_iw);
  }
  inline chi2_tau_t chi2_from_M2_PH(chi2_tau_t::const_view_type M2_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) {
    return chi2_from_M2<Chan_t::PH>(M2_tau, M_iw, G0_iw);
  }
  inline chi3_tau_t M3_conn_from_M3_PP(chi3_tau_t M3_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau,
                                       std::vector<matrix<M_tau_scalar_t>> M_hartree) {
    return M3_conn_from_M3<Chan_t::PP>(M3_tau, M_iw, G0_iw, M_tau, M_hartree);
  }
  inline chi3_tau_t M3_conn_from_M3_PH(chi3_tau_t M3_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau,
                                       std::vector<matrix<M_tau_scalar_t>> M_hartree) {
    return M3_conn_from_M3<Chan_t::PH>(M3_tau, M_iw, G0_iw, M_tau, M_hartree);
  }
  inline chi2_tau_t M2_from_M3_PP(chi3_tau_t M3_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau, std::vector<matrix<M_tau_scalar_t>> M_hartree,
                                  g_tau_cv_t G0_tau, int n_tau_M2) {
    return M2_from_M3<Chan_t::PP>(M3_tau, M_iw, G0_iw, M_tau, M_hartree, G0_tau, n_tau_M2);
  }
  inline chi2_tau_t M2_from_M3_PH(chi3_tau_t M3_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau, std::vector<matrix<M_tau_scalar_t>> M_hartree,
                                  g_tau_cv_t G0_tau, int n_tau_M2) {
    return M2_from_M3<Chan_t::PH>(M3_tau, M_iw, G0_iw, M_tau, M_hartree, G0_tau, n_tau_M2);
  }

} // namespace triqs_ctint
