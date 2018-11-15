#include "./M3pp_tau.hpp"

namespace triqs_ctint::measures {

  M3pp_tau::M3pp_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(std::move(G0_tau_)), tau_mesh{params_.beta, Fermion, params_.n_tau_M3} {

    // Construct Matsubara mesh
    gf_mesh<cartesian_product<imtime, imtime>> M3pp_tau_mesh{tau_mesh, tau_mesh};

    // Init measurement container and capture view
    results->M3pp_tau = make_block2_gf(M3pp_tau_mesh, params.gf_struct);
    M3pp_tau_.rebind(results->M3pp_tau.value());
    M3pp_tau_() = 0;
  }

  void M3pp_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Vectors containing the binned tau-values vec[bl][i]
    std::vector<std::vector<idx_t>> c_vec, cdag_vec;

    // The two tau meshes in one dimension
    auto const &G0_tau_mesh = G0_tau[0].mesh();
    auto const &M_tau_mesh  = tau_mesh;

    // Precompute binned tau-points
    for (auto &det : qmc_config.dets) {

      // Consider shifted time here as need in fourier transform (for the transform of unbarred index of M)
      auto x_to_mesh = [ beta = params.beta, &M_tau_mesh ](c_t const &c_i) { return idx_t{bin_to_mesh(beta - double(c_i.tau), M_tau_mesh), c_i.u}; };

      auto y_to_mesh = [ beta = params.beta, &G0_tau_mesh ](cdag_t const &cdag_j) {
        return idx_t{bin_to_mesh(beta - double(cdag_j.tau), G0_tau_mesh), cdag_j.u};
      };

      // Careful: Use the row and column indices of the matrix in their internal storage order
      c_vec.push_back(make_vector_from_range(transform(det.get_x_internal_order(), x_to_mesh)));
      cdag_vec.push_back(make_vector_from_range(transform(det.get_y_internal_order(), y_to_mesh)));
    }

    // The intermediate scattering matrix
    std::vector<matrix<dcomplex>> GM_vec(params.n_blocks()); // GM_vec[bl](u, i)

    // Calculate intermediate scattering matrix
    for (int bl : range(params.n_blocks())) {

      auto const &det = qmc_config.dets[bl];
      int det_size    = det.size();

      if (det.size() == 0) continue;

      auto const &cdag = cdag_vec[bl];
      int bl_size      = G0_tau[bl].target_shape()[0];
      auto G           = matrix<dcomplex>(bl_size, det_size);

      for (int b_u : range(bl_size))
        for (int j : range(det_size)) G(b_u, j) = G0_tau[bl][cdag[j].tau_idx](b_u, cdag[j].u);

      GM_vec[bl] = G * det.inverse_matrix_internal_order();
    }

    // Calculate M3pp
    for (int bl1 : range(params.n_blocks())) {

      int det1_size = qmc_config.dets[bl1].size();

      // Do not consider empty blocks
      if (det1_size == 0) continue;

      auto const &GM1 = GM_vec[bl1];
      auto const &c1  = c_vec[bl1];
      int bl1_size    = G0_tau[bl1].target_shape()[0];

      // Crossing term (equal blocks)
      auto &M3pp_tau = M3pp_tau_(bl1, bl1);
      // FIXME for( auto [k,l,i,j] : product_range(bl1_size, bl1_size, det1_size, det1_size) )
      for (int j : range(bl1_size))
        for (int l : range(bl1_size))
          for (int i : range(det1_size))
            for (int k : range(det1_size)) M3pp_tau[{c1[i].tau_idx, c1[k].tau_idx}](c1[i].u, j, c1[k].u, l) += -sign * GM1(l, i) * GM1(j, k);

      for (int bl2 : range(params.n_blocks())) {

        int det2_size = qmc_config.dets[bl2].size();

        // Do not consider empty blocks
        if (det2_size == 0) continue;

        auto const &GM2 = GM_vec[bl2];
        auto const &c2  = c_vec[bl2];
        int bl2_size    = G0_tau[bl2].target_shape()[0];
        auto &M3pp_tau  = M3pp_tau_(bl1, bl2);

        // Direct term
        for (int j : range(bl1_size))
          for (int l : range(bl2_size))
            for (int i : range(det1_size))
              for (int k : range(det2_size)) M3pp_tau[{c1[i].tau_idx, c2[k].tau_idx}](c1[i].u, j, c2[k].u, l) += sign * GM1(j, i) * GM2(l, k);
      }
    }
  }

  void M3pp_tau::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z           = mpi_all_reduce(Z, comm);
    M3pp_tau_   = mpi_all_reduce(M3pp_tau_, comm);
    double dtau = params.beta / (params.n_tau_M3 - 1);
    M3pp_tau_   = M3pp_tau_ / (Z * dtau * dtau);

    // Account for edge bins beeing smaller
    auto _ = all_t{};
    int n  = params.n_tau_M3 - 1;
    for (auto &M : M3pp_tau_) {
      M[0, _] *= 2.0;
      M[_, 0] *= 2.0;
      M[n, _] *= 2.0;
      M[_, n] *= 2.0;
    }
  }

} // namespace triqs_ctint::measures
