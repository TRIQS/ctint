// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "./qmc_config.hpp"
#include "./params.hpp"

namespace triqs_ctint {

  using triqs::operators::bilinear_type;

  /// Calculate the connected part of the two-particle Green function from M4_iw and M_iw
  chi4_iw_t G2_conn_from_M4(chi4_iw_t::const_view_type M4_iw, g_reg_iw_cv_t M_iw, g_reg_iw_cv_t G0_iw);
  /// Calculate the connected part of the two-particle Green function from M4pp_iw and M_iw
  chi4_iw_t G2pp_conn_from_M4pp(chi4_iw_t::const_view_type M4pp_iw, g_reg_iw_cv_t M_iw, g_reg_iw_cv_t G0_iw);
  /// Calculate the connected part of the two-particle Green function from M4pp_iw and M_iw
  chi4_iw_t G2ph_conn_from_M4ph(chi4_iw_t::const_view_type M4ph_iw, g_reg_iw_cv_t M_iw, g_reg_iw_cv_t G0_iw);

  /// Calculate the vertex function $F$ from G2_conn_iw and G_iw
  chi4_iw_t F_from_G2c(chi4_iw_t::const_view_type G2_conn_iw, g_reg_iw_cv_t G_iw);
  /// Calculate the vertex function $Fpp$ from G2pp_conn_iw and G_iw
  chi4_iw_t Fpp_from_G2pp_conn(chi4_iw_t::const_view_type G2pp_conn_iw, g_reg_iw_cv_t G_iw);
  /// Calculate the vertex function $Fph$ from G2ph_conn_iw and G_iw
  chi4_iw_t Fph_from_G2ph_conn(chi4_iw_t::const_view_type G2ph_conn_iw, g_reg_iw_cv_t G_iw);

  /// Calculate the two-particle Green function from G2_conn_iw and G_iw
  chi4_iw_t G2_from_G2c(chi4_iw_t::const_view_type G2_conn_iw, g_reg_iw_cv_t G_iw);
  /// Calculate the two-particle Green function from G2pp_conn_iw and G_iw
  chi4_iw_t G2pp_from_G2pp_conn(chi4_iw_t::const_view_type G2pp_conn_iw, g_reg_iw_cv_t G_iw);
  /// Calculate the two-particle Green function from G2ph_conn_iw and G_iw
  chi4_iw_t G2ph_from_G2ph_conn(chi4_iw_t::const_view_type G2ph_conn_iw, g_reg_iw_cv_t G_iw);

  /// Calculate the generalized ph susceptibility from G2ph_conn_iw and G_iw
  chi4_iw_t chi_tilde_ph_from_G2ph_conn(chi4_iw_t::const_view_type G2ph_conn_iw, g_reg_iw_cv_t G_iw);

  /// Calculate the $\chi_3$ function from the building blocks M3_iw and M_iw (DLR2D version)
  template <Chan_t Chan>
  chi3_dlr2d_iw_t chi3_from_M3(chi3_dlr2d_iw_cv_t M3_iw, g_reg_iw_cv_t M_iw, g_reg_iw_cv_t G0_iw, block_matrix_t const &dens_G,
                               block_matrix_t const &M_hartree) {

    double beta  = M_iw[0].mesh().beta();
    int n_blocks = M_iw.size();

    // Verify that input Green's function meshes cover all DLR2D frequencies
    long n_iw_G0 = static_cast<long>(G0_iw[0].mesh().size()) / 2;
    long max_n   = M3_iw(0, 0).mesh().max_n();
    if (max_n >= n_iw_G0)
      TRIQS_RUNTIME_ERROR << "chi3_from_M3: G0_iw mesh (n_iw=" << n_iw_G0
                          << ") does not cover all DLR2D mesh frequencies (max_n=" << max_n << ")";

    // Connected part of M3
    chi3_dlr2d_iw_t M3_iw_conn = M3_iw;

    // Connected part of chi3
    chi3_dlr2d_iw_t chi3_iw = M3_iw;
    chi3_iw()               = 0.;

    // Temporary quantities
    g_reg_iw_t GM   = G0_iw * M_iw;
    g_reg_iw_t MG   = M_iw * G0_iw;
    g_reg_iw_t GMG  = G0_iw * M_iw * G0_iw;
    g_reg_iw_t G_iw = G0_iw + G0_iw * M_iw * G0_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        // Capture block-sizes
        int bl1_size = M3_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M3_iw(bl1, bl2).target_shape()[2];

        if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

          // Loop over DLR2D mesh points
          for (auto mp : M3_iw(bl1, bl2).mesh()) {
            auto [iw1, iw2] = mp.value();
            for (int i : range(bl1_size))
              for (int j : range(bl1_size))
                for (int k : range(bl2_size))
                  for (int l : range(bl2_size)) {
                    M3_iw_conn(bl1, bl2)[mp](i, j, k, l) =
                       M3_iw(bl1, bl2)[mp](i, j, k, l) - GM[bl1][iw1](j, i) * GM[bl2][iw2](l, k) + kronecker(bl1, bl2) * GM[bl1][iw1](l, i) * GM[bl2][iw2](j, k);
                  }
          }

          for (auto mp : chi3_iw(bl1, bl2).mesh()) {
            auto [iw1, iw2] = mp.value();
            for (int i : range(bl1_size))
              for (int j : range(bl1_size))
                for (int k : range(bl2_size))
                  for (int l : range(bl2_size)) {
                    for (int m : range(bl1_size))
                      for (int n : range(bl2_size))
                        chi3_iw(bl1, bl2)[mp](i, j, k, l) += G0_iw[bl1][iw1](m, i) * G0_iw[bl2][iw2](n, k) * M3_iw_conn(bl1, bl2)[mp](m, j, n, l);
                    // Disconnected part
                    chi3_iw(bl1, bl2)[mp](i, j, k, l) +=
                       G_iw[bl1][iw1](j, i) * G_iw[bl2][iw2](l, k) - kronecker(bl1, bl2) * G_iw[bl1][iw1](l, i) * G_iw[bl2][iw2](j, k);
                  }
          }

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

          auto km_GMG = make_zero_tail(GMG, 3);
          for (auto [km_bl, M_hartree_bl] : zip(km_GMG, M_hartree)) km_bl(2, ellipsis()) = M_hartree_bl;
          auto tail_GMG = fit_hermitian_tail(GMG, km_GMG).first;
          auto dens_GMG = density(GMG, tail_GMG);

          // Loop over DLR2D mesh points
          for (auto mp : M3_iw(bl1, bl2).mesh()) {
            auto [iw1, iw2] = mp.value();
            for (int i : range(bl1_size))
              for (int j : range(bl1_size))
                for (int k : range(bl2_size))
                  for (int l : range(bl2_size)) {
                    M3_iw_conn(bl1, bl2)[mp](i, j, k, l) = M3_iw(bl1, bl2)[mp](i, j, k, l)
                       - beta * kronecker(iw1.n, iw2.n) * M_iw[bl1][iw1](j, i) * dens_GMG[bl2](l, k)
                       + kronecker(bl1, bl2) * GM[bl1][iw1](l, i) * MG[bl2][iw2](j, k);
                  }
          }

          for (auto mp : chi3_iw(bl1, bl2).mesh()) {
            auto [iw1, iw2] = mp.value();
            for (int i : range(bl1_size))
              for (int j : range(bl1_size))
                for (int k : range(bl2_size))
                  for (int l : range(bl2_size)) {
                    for (int m : range(bl1_size))
                      for (int n : range(bl1_size))
                        chi3_iw(bl1, bl2)[mp](i, j, k, l) += G0_iw[bl1][iw1](m, i) * G0_iw[bl1][iw2](j, n) * M3_iw_conn(bl1, bl2)[mp](m, n, k, l);
                    // Disconnected part
                    chi3_iw(bl1, bl2)[mp](i, j, k, l) +=
                       beta * kronecker(iw1.n, iw2.n) * G_iw[bl1][iw1](j, i) * dens_G[bl2](l, k) - kronecker(bl1, bl2) * G_iw[bl1][iw1](l, i) * G_iw[bl2][iw2](j, k);
                  }
          }
        }
      }

    return chi3_iw;
  }

  /// Calculate the $\chi_3$ function from the building blocks M3_iw and M_iw
  template <Chan_t Chan>
  chi3_iw_t chi3_from_M3(chi3_iw_cv_t M3_iw, g_reg_iw_cv_t M_iw, g_reg_iw_cv_t G0_iw, block_matrix_t const &dens_G, block_matrix_t const &M_hartree) {

    double beta  = M_iw[0].mesh().beta();
    int n_blocks = M_iw.size();

    // Connected part of M3
    chi3_iw_t M3_iw_conn = M3_iw;

    // Connected part of chi3
    chi3_iw_t chi3_iw = M3_iw;
    chi3_iw()         = 0.;

    // Temporary quantities
    g_reg_iw_t GM   = G0_iw * M_iw;
    g_reg_iw_t MG   = M_iw * G0_iw;
    g_reg_iw_t GMG  = G0_iw * M_iw * G0_iw;
    g_reg_iw_t G_iw = G0_iw + G0_iw * M_iw * G0_iw;

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        // Capture block-sizes
        int bl1_size = M3_iw(bl1, bl2).target_shape()[0];
        int bl2_size = M3_iw(bl1, bl2).target_shape()[2];

        if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

          M3_iw_conn(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                - GM[bl1](iw1_)(j_, i_) * GM[bl2](iw2_)(l_, k_) + kronecker(bl1, bl2) * GM[bl1](iw1_)(l_, i_) * GM[bl2](iw2_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl2_size))
              chi3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                    + G0_iw[bl1](iw1_)(m, i_) * G0_iw[bl2](iw2_)(n, k_) * M3_iw_conn(bl1, bl2)(iw1_, iw2_)(m, j_, n, l_);

          // Disconnected part
          chi3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                + G_iw[bl1](iw1_)(j_, i_) * G_iw[bl2](iw2_)(l_, k_) - kronecker(bl1, bl2) * G_iw[bl1](iw1_)(l_, i_) * G_iw[bl2](iw2_)(j_, k_);

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

          auto km_GMG = make_zero_tail(GMG, 3);
          for (auto [km_bl, M_hartree_bl] : zip(km_GMG, M_hartree)) km_bl(2, ellipsis()) = M_hartree_bl;
          auto tail_GMG = fit_hermitian_tail(GMG, km_GMG).first;
          auto dens_GMG = density(GMG, tail_GMG);

          M3_iw_conn(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                - beta * kronecker(iw1_, iw2_) * M_iw[bl1](iw1_)(j_, i_) * dens_GMG[bl2](l_, k_)
                + kronecker(bl1, bl2) * GM[bl1](iw1_)(l_, i_) * MG[bl2](iw2_)(j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl1_size))
              chi3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                    + G0_iw[bl1](iw1_)(m, i_) * G0_iw[bl1](iw2_)(j_, n) * M3_iw_conn(bl1, bl2)(iw1_, iw2_)(m, n, k_, l_);

          // Disconnected part
          chi3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                + beta * kronecker(iw1_, iw2_) * G_iw[bl1](iw1_)(j_, i_) * dens_G[bl2](l_, k_)
                - kronecker(bl1, bl2) * G_iw[bl1](iw1_)(l_, i_) * G_iw[bl2](iw2_)(j_, k_);
        } else if constexpr (Chan == Chan_t::XPH) { // ===== Particle-hole-cross channel

          auto km_GMG = make_zero_tail(GMG, 3);
          for (auto [km_bl, M_hartree_bl] : zip(km_GMG, M_hartree)) km_bl(2, ellipsis()) = M_hartree_bl;
          auto tail_GMG = fit_hermitian_tail(GMG, km_GMG).first;
          auto dens_GMG = density(GMG, tail_GMG);

          M3_iw_conn(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << M3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                - GM[bl1](iw1_)(j_, i_) * MG[bl2](iw2_)(l_, k_)
                + kronecker(bl1, bl2) * beta * kronecker(iw1_, iw2_) * M_iw[bl1](iw1_)(l_, i_) * dens_GMG[bl2](j_, k_);

          for (int m : range(bl1_size))
            for (int n : range(bl1_size))
              chi3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                    + G0_iw[bl1](iw1_)(m, i_) * G0_iw[bl2](iw2_)(l_, n) * M3_iw_conn(bl1, bl2)(iw1_, iw2_)(m, j_, k_, n);

          // Disconnected part
          chi3_iw(bl1, bl2)(iw1_, iw2_)(i_, j_, k_, l_) << chi3_iw(bl1, bl2)[iw1_, iw2_](i_, j_, k_, l_)
                + G_iw[bl1](iw1_)(j_, i_) * G_iw[bl2](iw2_)(l_, k_)
                - kronecker(bl1, bl2) * beta * kronecker(iw1_, iw2_) * G_iw[bl1](iw1_)(l_, i_) * dens_GMG[bl2](j_, k_);
        }
      }

    return chi3_iw;
  }

  // Calculate the $\chi_{AB}$ function from chi2_tau (DLR mesh)
  template <Chan_t Chan>
  gf<mesh::dlr_imtime, matrix_valued> chiAB_from_chi2(chi2_tau_cv_t chi2_tau, gf_struct_t const &gf_struct,
                                                       std::vector<many_body_operator> const &A_op_vec,
                                                       std::vector<many_body_operator> const &B_op_vec) {

    using op_term_t = std::tuple<dcomplex, bilinear_type, std::pair<int, int>, std::pair<int, int>>;
    std::vector<std::vector<op_term_t>> A_vec;
    std::vector<std::vector<op_term_t>> B_vec;

    for (auto A : A_op_vec) A_vec.emplace_back(get_terms(A, gf_struct));
    for (auto B : B_op_vec) B_vec.emplace_back(get_terms(B, gf_struct));

    auto chiAB_tau = gf<mesh::dlr_imtime, matrix_valued>{chi2_tau(0, 0).mesh(), make_shape(A_vec.size(), B_vec.size())};

    for (auto [j, B] : enumerate(B_vec))
      for (auto &[coef_B, type_B, bl_pair_B, idx_pair_B] : B) {

        if (type_B != bilinear_type::cdag_c)
          TRIQS_RUNTIME_ERROR << "chiAB_from_chi2 only supports c†c operators, got anomalous operator in B";

        auto [idx_cdag_B, idx_c_B] = idx_pair_B;
        auto [bl_cdag_B, bl_c_B]   = bl_pair_B;

        for (auto [i, A] : enumerate(A_vec))
          for (auto &[coef_A, type_A, bl_pair_A, idx_pair_A] : A) {

            if (type_A != bilinear_type::cdag_c)
              TRIQS_RUNTIME_ERROR << "chiAB_from_chi2 only supports c†c operators, got anomalous operator in A";

            auto [idx_cdag_A, idx_c_A] = idx_pair_A;
            auto [bl_cdag_A, bl_c_A]   = bl_pair_A;

            if ((bl_cdag_A != bl_c_A) || (bl_cdag_B != bl_c_B)) {
              TRIQS_RUNTIME_ERROR << "Monomials with unequal blocks not implemented for chiAB_from_chi2";
            }

            auto chiAB_tau_ij       = slice_target_to_scalar(chiAB_tau, i, j);
            auto chi2_tau_AABB_ijkl = slice_target_to_scalar(chi2_tau(bl_c_A, bl_c_B), idx_cdag_A, idx_c_A, idx_cdag_B, idx_c_B);
            chiAB_tau_ij() += coef_A * coef_B * chi2_tau_AABB_ijkl;
          }
      }

    return chiAB_tau;
  }

  // For wrapping purposes
  inline chi3_iw_t chi3_from_M3_PP(chi3_iw_cv_t M3_iw, g_reg_iw_cv_t M_iw, g_reg_iw_cv_t G0_iw, block_matrix_t const &dens_G,
                                   block_matrix_t const &M_hartree) {
    return chi3_from_M3<Chan_t::PP>(M3_iw, M_iw, G0_iw, dens_G, M_hartree);
  }
  inline chi3_iw_t chi3_from_M3_PH(chi3_iw_cv_t M3_iw, g_reg_iw_cv_t M_iw, g_reg_iw_cv_t G0_iw, block_matrix_t const &dens_G,
                                   block_matrix_t const &M_hartree) {
    return chi3_from_M3<Chan_t::PH>(M3_iw, M_iw, G0_iw, dens_G, M_hartree);
  }
  inline gf<mesh::dlr_imtime, matrix_valued> chiAB_from_chi2_PP(chi2_tau_cv_t chi2pp_tau, gf_struct_t const &gf_struct,
                                                                std::vector<many_body_operator> const &A_op_vec,
                                                                std::vector<many_body_operator> const &B_op_vec) {
    return chiAB_from_chi2<Chan_t::PP>(chi2pp_tau, gf_struct, A_op_vec, B_op_vec);
  }
  inline gf<mesh::dlr_imtime, matrix_valued> chiAB_from_chi2_PH(chi2_tau_cv_t chi2ph_tau, gf_struct_t const &gf_struct,
                                                                 std::vector<many_body_operator> const &A_op_vec,
                                                                 std::vector<many_body_operator> const &B_op_vec) {
    return chiAB_from_chi2<Chan_t::PH>(chi2ph_tau, gf_struct, A_op_vec, B_op_vec);
  }

} // namespace triqs_ctint
