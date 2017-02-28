#include "./M3_iw.hpp"

namespace triqs_ctint::measures {

  template <Chan_t Chan>
  M3_iw<Chan>::M3_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, block_gf<imtime, matrix_valued> const &G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(G0_tau_) {

    // Construct Matsubara mesh
    gf_mesh<imfreq> iw_mesh{params.beta, Fermion, params.n_iw_M3};
    gf_mesh<cartesian_product<imfreq, imfreq>> M3_iw_mesh{iw_mesh, iw_mesh};

    // Init measurement container
    std::vector<std::vector<block_type>> v_iw;
    for (auto const &bl1 : params.gf_struct) {
      std::vector<block_type> temp_vec;
      for (auto const &bl2 : params.gf_struct)
        temp_vec.emplace_back(M3_iw_mesh, make_shape(bl1.second.size(), bl1.second.size(), bl2.second.size(), bl2.second.size()));
      v_iw.emplace_back(std::move(temp_vec));
    }

    // Capture view to corresponding container_set member
    if (Chan == Chan_t::PP) { // FIXME c++17 if constexpr
      results->M3pp_iw = make_block2_gf(params.block_names(), params.block_names(), v_iw);
      M3_iw_.rebind(*results->M3pp_iw);
    } else if (Chan == Chan_t::PH) {
      results->M3ph_iw = make_block2_gf(params.block_names(), params.block_names(), v_iw);
      M3_iw_.rebind(*results->M3ph_iw);
    } else if (Chan == Chan_t::XPH) {
      results->M3xph_iw = make_block2_gf(params.block_names(), params.block_names(), v_iw);
      M3_iw_.rebind(*results->M3xph_iw);
    }
    M3_iw_() = 0;

    // Create nfft buffers
    for (int b1 = 0; b1 < params.n_blocks(); ++b1) {
      std::vector<array<triqs::utility::nfft_buf_t<2>, 4>> temp_vec;
      for (int b2 = 0; b2 < params.n_blocks(); ++b2) {
        auto init_target_func = [&](int i, int j, int k, int l) {
          return triqs::utility::nfft_buf_t<2>{slice_target_to_scalar(M3_iw_(b1, b2), i, j, k, l).data(), params.nfft_buf_size, params.beta};
        };
        temp_vec.emplace_back(M3_iw_(b1, b2).target_shape(), init_target_func);
      }
      buf_vecvec.emplace_back(std::move(temp_vec));
    }
  }

  template <Chan_t Chan> void M3_iw<Chan>::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    for (int b1 = 0; b1 < params.n_blocks(); ++b1)
      foreach (qmc_config.dets[b1], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv1) {
        for (int b2 = 0; b2 < params.n_blocks(); ++b2)
          foreach (qmc_config.dets[b2], [&](c_t const &c_k, cdag_t const &cdag_l, auto const &Ginv2) {

            int b1_size = M3_iw_(b1, b2).target_shape()[0];
            int b2_size = M3_iw_(b1, b2).target_shape()[2];

            double tau_i            = double(c_i.tau);
            double beta_minus_tau_j = params.beta - double(cdag_j.tau);
            double tau_k            = double(c_k.tau);
            double beta_minus_tau_l = params.beta - double(cdag_l.tau);

            mc_weight_t Ginv1_x_Ginv2_x_sign = Ginv1 * Ginv2 * sign;

            // Particle-particle channel
            if (Chan == Chan_t::PP) { // FIXME c++17 if constexpr
              for (int b_u : range(b1_size))
                for (int d_u : range(b2_size)) {
                  auto G0_bj = G0_tau[b1][closest_mesh_pt(beta_minus_tau_j)](b_u, cdag_j.u); // FIXME real/imag
                  auto G0_dl = G0_tau[b2][closest_mesh_pt(beta_minus_tau_l)](d_u, cdag_l.u);
                  buf_vecvec[b1][b2](c_i.u, b_u, c_k.u, d_u).push_back({tau_i, tau_k}, G0_bj * G0_dl * Ginv1_x_Ginv2_x_sign);

                  if (b1 == b2) { // Negative part with j and l swapped
                    auto G0_bl = G0_tau[b1][closest_mesh_pt(beta_minus_tau_l)](b_u, cdag_l.u);
                    auto G0_dj = G0_tau[b2][closest_mesh_pt(beta_minus_tau_j)](d_u, cdag_j.u);
                    buf_vecvec[b1][b2](c_i.u, b_u, c_k.u, d_u).push_back({tau_i, tau_k}, -G0_bl * G0_dj * Ginv1_x_Ginv2_x_sign);
                  }
                }
            } else if (Chan == Chan_t::PH) { // Particle-hole channel
              for (int cbar_u : range(b2_size))
                for (int d_u : range(b2_size)) {
                  auto G0_kc = G0_tau[b2][closest_mesh_pt(c_k.tau)](c_k.u, cbar_u);
                  auto G0_dl = G0_tau[b2][closest_mesh_pt(beta_minus_tau_l)](d_u, cdag_l.u);
                  buf_vecvec[b1][b2](c_i.u, cdag_j.u, cbar_u, d_u).push_back({tau_i, beta_minus_tau_j}, G0_kc * G0_dl * Ginv1_x_Ginv2_x_sign);

                  if (b1 == b2) { // Negative part with j and l swapped
                    auto G0_dj = G0_tau[b2][closest_mesh_pt(beta_minus_tau_j)](d_u, cdag_j.u);
                    buf_vecvec[b1][b2](c_i.u, cdag_l.u, cbar_u, d_u).push_back({tau_i, beta_minus_tau_l}, -G0_kc * G0_dj * Ginv1_x_Ginv2_x_sign);
                  }
                }
            } else if (Chan == Chan_t::XPH) { // Particle-hole transverse channel
              for (int b_u : range(b1_size))
                for (int cbar_u : range(b2_size)) {
                  auto G0_bj = G0_tau[b1][closest_mesh_pt(beta_minus_tau_j)](b_u, cdag_j.u);
                  auto G0_kc = G0_tau[b2][closest_mesh_pt(c_k.tau)](c_k.u, cbar_u);
                  buf_vecvec[b1][b2](c_i.u, b_u, cbar_u, cdag_l.u).push_back({tau_i, beta_minus_tau_l}, G0_bj * G0_kc * Ginv1_x_Ginv2_x_sign);

                  if (b1 == b2) { // Negative part with j and l swapped
                    auto G0_bl = G0_tau[b1][closest_mesh_pt(beta_minus_tau_l)](b_u, cdag_l.u);
                    buf_vecvec[b1][b2](c_i.u, b_u, cbar_u, cdag_j.u).push_back({tau_i, beta_minus_tau_j}, -G0_bl * G0_kc * Ginv1_x_Ginv2_x_sign);
                  }
                }
            }

          })
            ;
      })
        ;
  }

  template <Chan_t Chan> void M3_iw<Chan>::collect_results(triqs::mpi::communicator const &comm) {
    // Flush remaining points in nfft buffers
    for (auto &buf_vec : buf_vecvec)
      for (auto &buf_arr : buf_vec)
        for (auto &buf : buf_arr) buf.flush();

    // Collect results and normalize
    Z      = mpi_all_reduce(Z, comm);
    M3_iw_ = mpi_all_reduce(M3_iw_, comm);
    M3_iw_ = M3_iw_ / Z;
  }

} // namespace triqs_ctint::measures
