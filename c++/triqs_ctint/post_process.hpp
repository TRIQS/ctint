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

  /// Calculate the generalized ph susceptibility from G2c_iw and G_iw
  chi4_iw_t chi_tilde_ph_from_G2c(chi4_iw_t::const_view_type G2c_iw, g_iw_cv_t G_iw, gf_struct_t const &gf_struct);

  /// Calculate the $\chi_3$ function from the building blocks M3_iw and M_iw
  template <Chan_t Chan> chi3_iw_t chi3_from_M3(chi3_iw_cv_t M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) {

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
  template <Chan_t Chan>
  chi2_tau_t chi2_from_M2(chi2_tau_cv_t M2_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, std::vector<matrix<M_tau_scalar_t>> const &M_hartree) {

    double beta  = M_iw[0].domain().beta;
    int n_blocks = M_iw.size();

    // Chi2 and its connected part
    chi2_tau_t chi2_tau      = M2_tau;
    chi2_tau_t chi2_tau_conn = M2_tau;

    // Temporary quantities
    g_iw_t GMG_iw = G0_iw * M_iw * G0_iw;
    g_iw_t G_iw   = G0_iw + GMG_iw;

    auto tau_mesh = make_adjoint_mesh(G_iw[0].mesh());

    auto km = make_zero_tail(GMG_iw, 3);
    for (auto [km_bl, M_hartree_bl] : zip(km, M_hartree)) km_bl(2, ellipsis()) = M_hartree_bl;
    auto tail = fit_hermitian_tail(GMG_iw, km).first;

    auto km_G = make_zero_tail(G_iw, 2);
    for (auto &km_bl : km_G) matrix_view<dcomplex>{km_bl(1, ellipsis())} = 1.0;
    auto tail_G = fit_hermitian_tail(G_iw, km_G).first;

#ifdef GTAU_IS_COMPLEX
    auto GMG_tau  = make_gf_from_fourier(GMG_iw, tau_mesh, tail);
    g_tau_t G_tau = make_gf_from_fourier(G_iw, tau_mesh, tail_G);
#else
    auto GMG_tau  = real(make_gf_from_fourier(GMG_iw, tau_mesh, tail));
    g_tau_t G_tau = real(make_gf_from_fourier(G_iw, tau_mesh, tail_G));
    if (!is_gf_real_in_tau(G_iw, 1e-8)) std::cerr << "WARNING: Assuming real G_tau, but found Imag(G(tau)) > 1e-8. Casting to Real.\n";
#endif

    auto dens_GMG = density(GMG_iw, tail);

    for (int bl1 : range(n_blocks))
      for (int bl2 : range(n_blocks)) {

        if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

          chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_) << M2_tau(bl1, bl2)[t_](i_, j_, k_, l_)
                - GMG_tau[bl1](beta - t_)(j_, i_) * GMG_tau[bl2](beta - t_)(l_, k_)
                + kronecker(bl1, bl2) * GMG_tau[bl1](beta - t_)(l_, i_) * GMG_tau[bl2](beta - t_)(j_, k_);

          // Add Disconnected part
          chi2_tau(bl1, bl2)(t_)(i_, j_, k_, l_) << chi2_tau_conn(bl1, bl2)[t_](i_, j_, k_, l_)
                + G_tau[bl1](beta - t_)(j_, i_) * G_tau[bl2](beta - t_)(l_, k_)
                - kronecker(bl1, bl2) * G_tau[bl1](beta - t_)(l_, i_) * G_tau[bl2](beta - t_)(j_, k_);

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

          chi2_tau_conn(bl1, bl2)(t_)(i_, j_, k_, l_) << M2_tau(bl1, bl2)[t_](i_, j_, k_, l_) - dens_GMG[bl1](j_, i_) * dens_GMG[bl2](l_, k_)
                - kronecker(bl1, bl2) * GMG_tau[bl1](beta - t_)(l_, i_) * GMG_tau[bl2](t_)(j_, k_);

          // Add Disconnected part, FIXME used measured density?
          chi2_tau(bl1, bl2)(t_)(i_, j_, k_, l_) << chi2_tau_conn(bl1, bl2)[t_](i_, j_, k_, l_)
                + density(G_iw[bl1])(j_, i_) * density(G_iw[bl2])(l_, k_)
                + kronecker(bl1, bl2) * G_tau[bl1](beta - t_)(l_, i_) * G_tau[bl2](t_)(j_, k_); // Sign-change from G_tau shift
        }
      }

    return chi2_tau;
  }

  // Calculate the $\chi_2$ function from the building blocks M2_tau and M_iw
  template <Chan_t Chan>
  gf<imtime, matrix_valued> chiAB_from_chi2(chi2_tau_cv_t chi2_tau, gf_struct_t const &gf_struct, std::vector<many_body_operator> const &A_op_vec,
                                            std::vector<many_body_operator> const &B_op_vec) {

    using op_term_t = std::tuple<dcomplex, std::pair<int, int>, std::pair<int, int>>;
    std::vector<std::vector<op_term_t>> A_vec;
    std::vector<std::vector<op_term_t>> B_vec;

    for (auto A : A_op_vec) A_vec.emplace_back(get_terms(A, gf_struct));
    for (auto B : B_op_vec) B_vec.emplace_back(get_terms(B, gf_struct));

    auto chiAB_tau = gf<imtime, matrix_valued>{chi2_tau(0, 0).mesh(), make_shape(A_vec.size(), B_vec.size())};

    for (auto [j, B] : enumerate(B_vec))
      for (auto &[coef_B, bl_pair_B, idx_pair_B] : B) {

        auto [idx_cdag_B, idx_c_B] = idx_pair_B;
        auto [bl_cdag_B, bl_c_B]   = bl_pair_B;

        for (auto [i, A] : enumerate(A_vec))
          for (auto &[coef_A, bl_pair_A, idx_pair_A] : A) {

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

  /// Calculate the chi2_tau from M3_tau and M_iw
  template <Chan_t Chan>
  chi2_tau_t M2_from_M3(chi3_tau_cv_t M3_tau, chi2_tau_cv_t M3_delta, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau,
                        std::vector<matrix<M_tau_scalar_t>> const &M_hartree, g_tau_cv_t G0_tau, int n_tau_M2, gf_struct_t const &gf_struct) {

    double beta  = G0_tau[0].domain().beta;
    int n_blocks = G0_tau.size();

    auto const &tau_mesh_M3 = M3_tau(0, 0).mesh();
    double dtau_M3          = std::get<0>(tau_mesh_M3).delta();
    int n_tau_M3            = std::get<0>(tau_mesh_M3).size();

    auto const &tau_mesh_M3_del = M3_delta(0, 0).mesh();
    double dtau_M3_del          = tau_mesh_M3_del.delta();
    double n_tau_M3_del         = tau_mesh_M3_del.size();

    auto tau_mesh_M2 = gf_mesh<imtime>{beta, Boson, n_tau_M2};
    auto M2_tau      = make_block2_gf(tau_mesh_M2, make_const_view(M3_tau));
    M2_tau()         = 0.0;

    // Set up MPI Parellelization
    // The following operation is memory bound,
    // so we perform the computation only on every n-th thread
    const int n_skip    = 2;
    auto comm           = triqs::mpi::communicator();
    const int mpi_rank  = comm.rank();
    const int thread_id = mpi_rank / n_skip;
    const int n_threads = 1 + (comm.size() - 1) / n_skip;

    // Some Helper variables
    double weight_del = sqrt(dtau_M3_del);
    double sqrt_2     = sqrt(0.5);

    for (auto &t : tau_mesh_M2) {

      // MPI Parellelization
      if (mpi_rank % n_skip != 0) continue;
      int job_id = t.linear_index();
      if (job_id % n_threads != thread_id) continue;

      // ==== Precalculate all relevant G0_tau interpolation points

      // Preallocate the G0 vectors
      auto G0_ti_minus_t     = std::vector<array<dcomplex, 3>>{};
      auto G0_t_minus_ti     = std::vector<array<dcomplex, 3>>{};
      auto G0_ti_minus_t_del = std::vector<array<dcomplex, 3>>{};
      auto G0_t_minus_ti_del = std::vector<array<dcomplex, 3>>{};

      for (int bl : range(n_blocks)) {

        auto mat_shape = G0_tau[bl].target_shape();
        G0_ti_minus_t.emplace_back(mat_shape[0], mat_shape[1], n_tau_M3);
        G0_t_minus_ti.emplace_back(mat_shape[0], mat_shape[1], n_tau_M3);
        G0_ti_minus_t_del.emplace_back(mat_shape[0], mat_shape[1], n_tau_M3_del);
        G0_t_minus_ti_del.emplace_back(mat_shape[0], mat_shape[1], n_tau_M3_del);

        for (auto const &ti : std::get<0>(tau_mesh_M3)) {

          auto [s1, d_ti_t] = cyclic_difference(ti, t);
          auto [s2, d_t_ti] = cyclic_difference(t, ti);

          // Treat certain equal-time cases seperately
          if (ti.index() == n_tau_M3 - 1 && t.index() == n_tau_M2 - 1) {
            s1     = -1.0;
            d_ti_t = beta - 1e-14;
          }
          if (ti.index() == 0 && t.index() == 0) {
            s2     = -1.0;
            d_t_ti = beta - 1e-14;
          }

          //Account for M3 edge bins beeing smaller
          if (ti.linear_index() == 0 or ti.linear_index() == n_tau_M3 - 1) {
            s1 *= 0.5;
            s2 *= 0.5;
          }

          G0_ti_minus_t[bl](range(), range(), ti.linear_index()) = dtau_M3 * s1 * G0_tau[bl](d_ti_t);
          G0_t_minus_ti[bl](range(), range(), ti.linear_index()) = dtau_M3 * s2 * G0_tau[bl](d_t_ti);
        }

        // Now for the Mesh of M3_del
        for (auto const &ti : tau_mesh_M3_del) {

          auto [s1, d_ti_t] = cyclic_difference(ti, t);
          auto [s2, d_t_ti] = cyclic_difference(t, ti);

          // Carefully treat certain equal-time cases
          if (ti.index() == n_tau_M3_del - 1 && t.index() == n_tau_M2 - 1) {
            s1     = -1.0;
            d_ti_t = beta - 1e-14;
          }
          if (ti.index() == 0 && t.index() == 0) {
            s2     = -1.0;
            d_t_ti = beta - 1e-14;
          }

          // Account for M3 delta edge bins beeing smaller
          if (ti.linear_index() == 0 or ti.linear_index() == n_tau_M3_del) {
            s1 *= sqrt_2;
            s2 *= sqrt_2;
          }

          G0_ti_minus_t_del[bl](range(), range(), ti.linear_index()) = weight_del * s1 * G0_tau[bl](d_ti_t);
          G0_t_minus_ti_del[bl](range(), range(), ti.linear_index()) = weight_del * s2 * G0_tau[bl](d_t_ti);
        }
      }
      // ==== End Precalculation

      for (auto [bl1, bl2] : product_range(n_blocks, n_blocks)) {

        // Capture block-sizes
        int bl1_size = M3_tau(bl1, bl2).target_shape()[0];
        int bl2_size = M3_tau(bl1, bl2).target_shape()[2];

        auto M2 = M2_tau(bl1, bl2)[t];

        if constexpr (Chan == Chan_t::PP) { // =====  Particle-particle channel

          auto arr_GG = array<dcomplex, 5>(bl1_size, bl1_size, bl2_size, bl2_size, n_tau_M3_del);
          for (auto [m, i, n, k] : product_range(bl1_size, bl1_size, bl2_size, bl2_size)) {
            arr_GG(m, i, n, k, range()) = G0_ti_minus_t_del[bl1](m, i, range()) * G0_ti_minus_t_del[bl2](n, k, range());
          }

          for (auto [m, j, n, l] : product_range(bl1_size, bl1_size, bl2_size, bl2_size)) {

            // We have to make a copy so that M3 is contiguous in memory
            auto M3_mjnl = matrix<dcomplex>{slice_target_to_scalar(M3_tau(bl1, bl2), m, j, n, l).data()};
            for (auto [i, k] : product_range(bl1_size, bl2_size)) {
              auto G1_mi = vector_view<dcomplex>(G0_ti_minus_t[bl1](m, i, range()));
              auto G2_nk = vector_view<dcomplex>(G0_ti_minus_t[bl2](n, k, range()));
              M2(i, j, k, l) += dot(G1_mi, M3_mjnl * G2_nk);
            }

            // We treat the delta-contribution seperately
            auto M3_del_mjnl = vector<dcomplex>(slice_target_to_scalar(M3_delta(bl1, bl2), m, j, n, l).data());
            for (auto [i, k] : product_range(bl1_size, bl2_size)) {
              M2(i, j, k, l) += dot(M3_del_mjnl, vector_view<dcomplex>(arr_GG(m, i, n, k, range())));
            }
          }

        } else if constexpr (Chan == Chan_t::PH) { // ===== Particle-hole channel

          auto arr_GG = array<dcomplex, 5>(bl1_size, bl1_size, bl1_size, bl1_size, n_tau_M3_del);
          for (auto [m, i, j, n] : product_range(bl1_size, bl1_size, bl1_size, bl1_size)) {
            arr_GG(m, i, j, n, range()) = G0_ti_minus_t_del[bl1](m, i, range()) * G0_t_minus_ti_del[bl1](j, n, range());
          }

          for (auto [m, n, k, l] : product_range(bl1_size, bl1_size, bl2_size, bl2_size)) {

            // We have to make a copy so that M3 is contiguous in memory
            auto M3_mnkl = matrix<dcomplex>{slice_target_to_scalar(M3_tau(bl1, bl2), m, n, k, l).data()};
            for (auto [i, j] : product_range(bl1_size, bl1_size)) {
              auto G1_mi = vector_view<dcomplex>(G0_ti_minus_t[bl1](m, i, range()));
              auto G2_jn = vector_view<dcomplex>(G0_t_minus_ti[bl1](j, n, range()));
              M2(i, j, k, l) += dot(G1_mi, M3_mnkl * G2_jn);
            }

            // We treat the delta-contribution seperately
            auto M3_del_mnkl = vector<dcomplex>(slice_target_to_scalar(M3_delta(bl1, bl2), m, n, k, l).data());
            for (auto [i, j] : product_range(bl1_size, bl1_size)) {
              M2(i, j, k, l) += dot(M3_del_mnkl, vector_view<dcomplex>(arr_GG(m, i, j, n, range())));
            }
          }
        }
      }
    }

    M2_tau() = mpi_all_reduce(M2_tau, comm);
    return M2_tau;
  }

  // For wrapping purposes
  inline chi3_iw_t chi3_from_M3_PP(chi3_iw_cv_t M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) { return chi3_from_M3<Chan_t::PP>(M3_iw, M_iw, G0_iw); }
  inline chi3_iw_t chi3_from_M3_PH(chi3_iw_cv_t M3_iw, g_iw_cv_t M_iw, g_iw_cv_t G0_iw) { return chi3_from_M3<Chan_t::PH>(M3_iw, M_iw, G0_iw); }
  inline chi2_tau_t chi2_from_M2_PP(chi2_tau_cv_t M2_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, std::vector<matrix<M_tau_scalar_t>> const &M_hartree) {
    return chi2_from_M2<Chan_t::PP>(M2_tau, M_iw, G0_iw, M_hartree);
  }
  inline chi2_tau_t chi2_from_M2_PH(chi2_tau_cv_t M2_tau, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, std::vector<matrix<M_tau_scalar_t>> const &M_hartree) {
    return chi2_from_M2<Chan_t::PH>(M2_tau, M_iw, G0_iw, M_hartree);
  }
  inline gf<imtime, matrix_valued> chiAB_from_chi2_PP(chi2_tau_cv_t chi2pp_tau, gf_struct_t const &gf_struct,
                                                      std::vector<many_body_operator> const &A_op_vec,
                                                      std::vector<many_body_operator> const &B_op_vec) {
    return chiAB_from_chi2<Chan_t::PP>(chi2pp_tau, gf_struct, A_op_vec, B_op_vec);
  }
  inline gf<imtime, matrix_valued> chiAB_from_chi2_PH(chi2_tau_cv_t chi2ph_tau, gf_struct_t const &gf_struct,
                                                      std::vector<many_body_operator> const &A_op_vec,
                                                      std::vector<many_body_operator> const &B_op_vec) {
    return chiAB_from_chi2<Chan_t::PH>(chi2ph_tau, gf_struct, A_op_vec, B_op_vec);
  }
  inline chi2_tau_t M2_from_M3_PP(chi3_tau_t M3pp_tau, chi2_tau_t M3pp_delta, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau,
                                  std::vector<matrix<M_tau_scalar_t>> const &M_hartree, g_tau_cv_t G0_tau, int n_tau_M2,
                                  gf_struct_t const &gf_struct) {
    return M2_from_M3<Chan_t::PP>(M3pp_tau, M3pp_delta, M_iw, G0_iw, M_tau, M_hartree, G0_tau, n_tau_M2, gf_struct);
  }
  inline chi2_tau_t M2_from_M3_PH(chi3_tau_t M3ph_tau, chi2_tau_t M3ph_delta, g_iw_cv_t M_iw, g_iw_cv_t G0_iw, g_tau_cv_t M_tau,
                                  std::vector<matrix<M_tau_scalar_t>> const &M_hartree, g_tau_cv_t G0_tau, int n_tau_M2,
                                  gf_struct_t const &gf_struct) {
    return M2_from_M3<Chan_t::PH>(M3ph_tau, M3ph_delta, M_iw, G0_iw, M_tau, M_hartree, G0_tau, n_tau_M2, gf_struct);
  }

} // namespace triqs_ctint
