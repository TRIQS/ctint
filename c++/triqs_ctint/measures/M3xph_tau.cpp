#include "./M3xph_tau.hpp"

namespace triqs_ctint::measures {

  M3xph_tau::M3xph_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(std::move(G0_tau_)), tau_mesh{params_.beta, Fermion, params_.n_tau_M3} {

    // Construct Matsubara mesh
    mesh::prod<imtime, imtime> M3xph_tau_mesh{tau_mesh, tau_mesh};

    // Init measurement container for M3xph and capture view
    results->M3xph_tau = make_block2_gf(M3xph_tau_mesh, params.gf_struct);
    M3xph_tau_.rebind(results->M3xph_tau.value());
    M3xph_tau_() = 0;

    // Init measurement container for equal-time component of M3xph
    results->M3xph_delta = make_block2_gf(G0_tau[0].mesh(), params.gf_struct);
    M3xph_delta_.rebind(*results->M3xph_delta);
    M3xph_delta_() = 0;
  }

  void M3xph_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Vectors containing the binned tau-values vec[bl][i] (binned to G0 and M3 mesh)
    std::vector<std::vector<idx_t>> c_vec_G0, cdag_vec_G0;
    std::vector<std::vector<idx_t>> c_vec_M, cdag_vec_M;

    // The two tau meshes in one dimension
    auto const &G0_tau_mesh = G0_tau[0].mesh();
    auto const &M_tau_mesh  = tau_mesh;

    // Precompute binned tau-points for different meshes
    for (auto &det : qmc_config.dets) {

      auto x_to_G0_mesh = [&G0_tau_mesh](c_t const &c_i) { return idx_t{G0_tau_mesh.to_idx(closest_mesh_pt(double(c_i.tau))), c_i.u, c_i.tau}; };
      auto y_to_G0_mesh = [beta = params.beta, &G0_tau_mesh](cdag_t const &cdag_j) {
        return idx_t{G0_tau_mesh.to_idx(closest_mesh_pt(beta - double(cdag_j.tau))), cdag_j.u, cdag_j.tau};
      };

      auto x_to_M_mesh = [&M_tau_mesh](c_t const &c_i) { return idx_t{M_tau_mesh.to_idx(closest_mesh_pt(double(c_i.tau))), c_i.u, c_i.tau}; };
      auto y_to_M_mesh = [&M_tau_mesh](cdag_t const &cdag_j) {
        return idx_t{M_tau_mesh.to_idx(closest_mesh_pt(double(cdag_j.tau))), cdag_j.u, cdag_j.tau};
      };

      // Careful: x and y vectors have to be used in internal storage order
      c_vec_G0.push_back(make_vector_from_range(transform(det.get_x_internal_order(), x_to_G0_mesh)));
      cdag_vec_G0.push_back(make_vector_from_range(transform(det.get_y_internal_order(), y_to_G0_mesh)));
      c_vec_M.push_back(make_vector_from_range(transform(det.get_x_internal_order(), x_to_M_mesh)));
      cdag_vec_M.push_back(make_vector_from_range(transform(det.get_y_internal_order(), y_to_M_mesh)));
    }

    // The intermediate scattering matrices
    std::vector<matrix<dcomplex>> M_vec(params.n_blocks());   // M[bl](j, i) FIXME matrix_view does not work!
    std::vector<matrix<dcomplex>> GM_vec(params.n_blocks());  // GM[bl](u, i)
    std::vector<matrix<dcomplex>> MG_vec(params.n_blocks());  // MG[bl](j, u)
    std::vector<matrix<dcomplex>> GMG_vec(params.n_blocks()); // GMG[bl](u, u)

    // Calculate intermediate scattering matrices
    for (int bl : range(params.n_blocks())) {

      auto const &det = qmc_config.dets[bl];
      int det_size    = det.size();

      if (det.size() == 0) continue;

      auto const &c    = c_vec_G0[bl];
      auto const &cdag = cdag_vec_G0[bl];
      int bl_size      = G0_tau[bl].target_shape()[0];
      auto G_left      = matrix<dcomplex>(bl_size, det_size);
      auto G_right     = matrix<dcomplex>(det_size, bl_size);

      for (int u : range(bl_size))
        for (int i : range(det_size)) {
          G_left(u, i)  = -G0_tau[bl][cdag[i].tau_idx](u, cdag[i].u);
          G_right(i, u) = G0_tau[bl][c[i].tau_idx](c[i].u, u);
        }
      M_vec[bl]   = det.inverse_matrix_internal_order();
      GM_vec[bl]  = G_left * M_vec[bl];
      MG_vec[bl]  = M_vec[bl] * G_right;
      GMG_vec[bl] = GM_vec[bl] * G_right;
    }

    // Calculate M3xph
    for (int bl1 : range(params.n_blocks())) {

      int det1_size = qmc_config.dets[bl1].size();

      // Do not consider empty blocks
      if (det1_size == 0) continue;

      auto const &M     = M_vec[bl1];
      auto const &GM    = GM_vec[bl1];
      auto const &GMG   = GMG_vec[bl1];
      int bl1_size      = G0_tau[bl1].target_shape()[0];
      auto const &c1    = c_vec_M[bl1];
      auto const &cdag1 = cdag_vec_M[bl1];

      // Crossing term (equal blocks)
      auto &M3xph_tau   = M3xph_tau_(bl1, bl1);
      auto &M3xph_delta = M3xph_delta_(bl1, bl1);

      for (auto [i, j, k, l] : product_range(det1_size, bl1_size, bl1_size, det1_size)) {
        // Take care of equal-time peak separately
        if (c1[i].tau_pt == cdag1[l].tau_pt) {
          M3xph_delta[c_vec_G0[bl1][i].tau_idx](c1[i].u, j, k, cdag1[l].u) += -sign * M(l, i) * GMG(j, k);
        } else {
          // Since the crossing term is negative by itself, we get a negative sign here
          M3xph_tau[c1[i].tau_idx, cdag1[l].tau_idx](c1[i].u, j, k, cdag1[l].u) += -sign * M(l, i) * GMG(j, k);
        }
      }

      for (int bl2 : range(params.n_blocks())) {

        int det2_size = qmc_config.dets[bl2].size();

        // Do not consider empty blocks
        if (det2_size == 0) continue;

        auto const &MG    = MG_vec[bl2];
        int bl2_size      = G0_tau[bl2].target_shape()[0];
        auto const &cdag2 = cdag_vec_M[bl2];
        auto &M3xph_tau   = M3xph_tau_(bl1, bl2);
        auto &M3xph_delta = M3xph_delta_(bl1, bl2);

        // Direct term
        for (auto [i, j, k, l] : product_range(det1_size, bl1_size, bl2_size, det2_size)) {
          // Take care of equal-time peak separately
          if (c1[i].tau_pt == cdag2[l].tau_pt) {
            M3xph_delta[c_vec_G0[bl1][i].tau_idx](c1[i].u, j, k, cdag2[l].u) += sign * GM(j, i) * MG(l, k);
          } else {
            M3xph_tau[c1[i].tau_idx, cdag2[l].tau_idx](c1[i].u, j, k, cdag2[l].u) += sign * GM(j, i) * MG(l, k);
          }
        }
      }
    }
  }

  void M3xph_tau::collect_results(mpi::communicator const &comm) {
    // Collect results
    Z            = mpi::all_reduce(Z, comm);
    M3xph_tau_   = mpi::all_reduce(M3xph_tau_, comm);
    M3xph_delta_ = mpi::all_reduce(M3xph_delta_, comm);

    // Normalize
    int n       = params.n_tau_M3 - 1;
    double dtau = params.beta / n;
    M3xph_tau_  = M3xph_tau_ / (Z * dtau * dtau);

    int n_del       = params.n_tau - 1;
    double dtau_del = params.beta / n_del;
    M3xph_delta_    = M3xph_delta_ / (Z * dtau_del);

    // Account for edge bins beeing smaller
    auto _ = all_t{};
    for (auto [M, M_del] : zip(M3xph_tau_, M3xph_delta_)) {
      M[0, _] *= 2.0;
      M[_, 0] *= 2.0;
      M[n, _] *= 2.0;
      M[_, n] *= 2.0;
      M_del[0] *= 2.0;
      M_del[n_del] *= 2.0;
    }
  }

} // namespace triqs_ctint::measures
