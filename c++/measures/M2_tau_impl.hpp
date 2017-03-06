#include "./M2_tau.hpp"

namespace triqs_ctint::measures {

  template <Chan_t Chan>
  M2_tau<Chan>::M2_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results,
                       block_gf<imtime, matrix_valued> const &G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(G0_tau_) {

    // Construct Matsubara mesh
    gf_mesh<imtime> tau_mesh{params.beta, Boson, params.n_tau_M2};

    // Init measurement container and capture view // FIXME c++17 if constexpr
    if (Chan == Chan_t::PP) {
      results->M2pp_tau = make_block2_gf(tau_mesh, params.gf_struct);
      M2_tau_.rebind(*results->M2pp_tau);
    } else if (Chan == Chan_t::PH) {
      results->M2ph_tau = make_block2_gf(tau_mesh, params.gf_struct);
      M2_tau_.rebind(*results->M2ph_tau);
    } else if (Chan == Chan_t::XPH) {
      results->M2xph_tau = make_block2_gf(tau_mesh, params.gf_struct);
      M2_tau_.rebind(*results->M2xph_tau);
    }
    M2_tau_() = 0;
  }

  template <Chan_t Chan> void M2_tau<Chan>::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    for (int b1 : range(params.n_blocks()))
      foreach (qmc_config.dets[b1], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv1) {
        for (int b2 : range(params.n_blocks()))
          foreach (qmc_config.dets[b2], [&](c_t const &c_k, cdag_t const &cdag_l, auto const &Ginv2) {

            int b1_size = M2_tau_(b1, b2).target_shape()[0]; // FIXME Change gf_struct ...
            int b2_size = M2_tau_(b1, b2).target_shape()[2];

            double tau_i = double(c_i.tau);
            double tau_j = double(cdag_j.tau);
            double tau_k = double(c_k.tau);
            double tau_l = double(cdag_l.tau);

            auto beta                 = params.beta;
            auto Ginv1_x_Ginv2_x_sign = Ginv1 * Ginv2 * sign;

            auto G0 =
               [this](int b, double tau_c, double tau_cdag, int u_c, int u_cdag) {
                 if (tau_c < tau_cdag) return -G0_tau[b][closest_mesh_pt(tau_cdag - tau_c)](u_c, u_cdag);
                 return G0_tau[b][closest_mesh_pt(tau_c - tau_cdag)](u_c, u_cdag);
               }; 

            // Particle-particle channel
            for (auto &tau_pt : M2_tau_(b1,b2).mesh()) {
              double tau = double(tau_pt);

              for (int abar_u : range(b1_size))
                for (int b_u : range(b1_size))
                  for (int cbar_u : range(b2_size))
                    for (int d_u : range(b2_size))
                      if (Chan == Chan_t::PP) { // Particle-particle channel //FIXME c++17 if constexpr
                        auto G0_ia = G0(b1, tau_i, tau, c_i.u, abar_u);
                        auto G0_kc = G0(b2, tau_k, tau, c_k.u, cbar_u);
                        auto G0_bj = G0_tau[b1][closest_mesh_pt(beta - tau_j)](b_u, cdag_j.u);
                        auto G0_dl = G0_tau[b2][closest_mesh_pt(beta - tau_l)](d_u, cdag_l.u);
                        M2_tau_(b1, b2)[tau_pt](abar_u, b_u, cbar_u, d_u) += G0_ia * G0_kc * G0_bj * G0_dl * Ginv1_x_Ginv2_x_sign;
                        if (b1 == b2) {
                          auto G0_bl = G0_tau[b1][closest_mesh_pt(beta - tau_l)](b_u, cdag_l.u);
                          auto G0_dj = G0_tau[b2][closest_mesh_pt(beta - tau_j)](d_u, cdag_j.u);
                          M2_tau_(b1, b2)[tau_pt](abar_u, b_u, cbar_u, d_u) += -G0_ia * G0_kc * G0_bl * G0_dj * Ginv1_x_Ginv2_x_sign;
                        }
                      } else if (Chan == Chan_t::PH) { // Particle-hole channel
                        auto G0_ia = G0(b1, tau_i, tau, c_i.u, abar_u);
                        auto G0_kc = G0_tau[b2][closest_mesh_pt(c_k.tau)](c_k.u, cbar_u);
                        auto G0_bj = G0(b1, tau, tau_j, b_u, cdag_j.u);
                        auto G0_dl = G0_tau[b2][closest_mesh_pt(beta - tau_l)](d_u, cdag_l.u);
                        M2_tau_(b1, b2)[tau_pt](abar_u, b_u, cbar_u, d_u) += -G0_ia * G0_kc * G0_bj * G0_dl * Ginv1_x_Ginv2_x_sign;
                        if (b1 == b2) {
                          auto G0_bl = G0(b1, tau, tau_l, b_u, cdag_l.u);
                          auto G0_dj = G0_tau[b2][closest_mesh_pt(beta - tau_j)](d_u, cdag_j.u);
                          M2_tau_(b1, b2)[tau_pt](abar_u, b_u, cbar_u, d_u) += G0_ia * G0_kc * G0_bl * G0_dj * Ginv1_x_Ginv2_x_sign;
                        }
                      } else if (Chan == Chan_t::XPH) { // Particle-hole transverse channel
                        auto G0_ia = G0(b1, tau_i, tau, c_i.u, abar_u);
                        auto G0_kc = G0_tau[b2][closest_mesh_pt(c_k.tau)](c_k.u, cbar_u);
                        auto G0_bj = G0_tau[b1][closest_mesh_pt(beta - tau_j)](b_u, cdag_j.u);
                        auto G0_dl = G0(b2, tau, tau_l, d_u, cdag_l.u);
                        M2_tau_(b1, b2)[tau_pt](abar_u, b_u, cbar_u, d_u) += -G0_ia * G0_kc * G0_bj * G0_dl * Ginv1_x_Ginv2_x_sign;
                        if (b1 == b2) {
                          auto G0_bl = G0_tau[b1][closest_mesh_pt(beta - tau_l)](b_u, cdag_l.u);
                          auto G0_dj = G0(b2, tau, tau_j, d_u, cdag_j.u);
                          M2_tau_(b1, b2)[tau_pt](abar_u, b_u, cbar_u, d_u) += G0_ia * G0_kc * G0_bl * G0_dj * Ginv1_x_Ginv2_x_sign;
                        }
                      }
            }

          })
            ;
      })
        ;
  }

  template <Chan_t Chan> void M2_tau<Chan>::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z       = mpi_all_reduce(Z, comm);
    M2_tau_ = mpi_all_reduce(M2_tau_, comm);
    M2_tau_ = M2_tau_ / Z;
  }

} // namespace triqs_ctint::measures
