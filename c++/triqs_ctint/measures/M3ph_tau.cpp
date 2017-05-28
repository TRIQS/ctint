#include "./M3ph_tau.hpp"
#include "../itertools.hpp"

namespace triqs_ctint::measures {

  M3ph_tau::M3ph_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(G0_tau_), tau_mesh{params_.beta, Fermion, params_.n_tau_M3} {

    // Construct Matsubara mesh
    gf_mesh<cartesian_product<imtime, imtime>> M3ph_tau_mesh{tau_mesh, tau_mesh};

    // Init measurement container and capture view
    results->M3ph_tau = make_block2_gf(M3ph_tau_mesh, params.gf_struct);
    M3ph_tau_.rebind(*results->M3ph_tau);
    M3ph_tau_() = 0;
  }

  void M3ph_tau::accumulate(mc_weight_t sign) {
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

      auto x_to_G0_mesh = [&G0_tau_mesh](c_t const &c_i) { return idx_t{bin_to_mesh(double(c_i.tau), G0_tau_mesh), c_i.u}; };
      auto y_to_G0_mesh = [ beta = params.beta, &G0_tau_mesh ](cdag_t const &cdag_j) {
        return idx_t{bin_to_mesh(beta - double(cdag_j.tau), G0_tau_mesh), cdag_j.u};
      };

      // Careful: the unbarred (x) index of M has to be transformed with a negative time in the Fourier transform later -> shift
      auto x_to_M_mesh = [ beta = params.beta, &M_tau_mesh ](c_t const &c_i) {
        return idx_t{bin_to_mesh(beta - double(c_i.tau), M_tau_mesh), c_i.u};
      };
      auto y_to_M_mesh = [&M_tau_mesh](cdag_t const &cdag_j) { return idx_t{bin_to_mesh(double(cdag_j.tau), M_tau_mesh), cdag_j.u}; };

      // Careful: x and y vectors have to be used in internal storage order
      c_vec_G0.push_back(make_vector_from_range(transform(det.get_x_internal_order(), x_to_G0_mesh)));
      cdag_vec_G0.push_back(make_vector_from_range(transform(det.get_y_internal_order(), y_to_G0_mesh)));
      c_vec_M.push_back(make_vector_from_range(transform(det.get_x_internal_order(), x_to_M_mesh)));
      cdag_vec_M.push_back(make_vector_from_range(transform(det.get_y_internal_order(), y_to_M_mesh)));
    }

    // The intermediate scattering matrices
    std::vector<matrix<dcomplex>> M_vec;   // M[bl](j, i) FIXME matrix_view does not work!
    std::vector<matrix<dcomplex>> GM_vec;  // GM[bl](u, i)
    std::vector<matrix<dcomplex>> MG_vec;  // MG[bl](j, u)
    std::vector<matrix<dcomplex>> GMG_vec; // GMG[bl](u, u)

    // Calculate intermediate scattering matrices
    for (int bl : range(params.n_blocks())) {

      auto const &c    = c_vec_G0[bl];
      auto const &cdag = cdag_vec_G0[bl];
      auto const &det  = qmc_config.dets[bl];
      int bl_size      = G0_tau[bl].target_shape()[0];
      int det_size     = det.size();
      auto G_left      = matrix<dcomplex>(bl_size, det_size);
      auto G_right     = matrix<dcomplex>(det_size, bl_size);

      for (int u : range(bl_size))
        for (int i : range(det_size)) {
          G_left(u, i)  = -G0_tau[bl][cdag[i].tau_idx](u, cdag[i].u);
          G_right(i, u) = G0_tau[bl][c[i].tau_idx](c[i].u, u);
        }
      M_vec.emplace_back(det.inverse_matrix_internal_order());
      GM_vec.push_back(G_left * M_vec[bl]);
      MG_vec.push_back(M_vec[bl] * G_right);
      GMG_vec.push_back(GM_vec[bl] * G_right);
    }

    // Calculate M3ph
    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        auto const &M     = M_vec[bl1];
        auto const &GMG   = GMG_vec[bl2];
        auto const &GM    = GM_vec[bl1];
        auto const &MG    = MG_vec[bl2];
        auto const &c1    = c_vec_M[bl1];
        auto const &cdag1 = cdag_vec_M[bl1];
        int det1_size     = qmc_config.dets[bl1].size();
        int bl2_size      = G0_tau[bl2].target_shape()[0];
        auto &M3ph_tau    = M3ph_tau_(bl1, bl2);

        for (int k : range(bl2_size))
          for (int l : range(bl2_size))
            for (int i : range(det1_size))
              for (int j : range(det1_size))
                // We have to add an overall minus sign to account for the fourier transform of the right M index
                // FIXME Adjust multidimensional fourier transform to allow for sperarate e^{-iwt} transforms
                M3ph_tau[{c1[i].tau_idx, cdag1[j].tau_idx}](c1[i].u, cdag1[j].u, k, l) +=
                   -sign * (M(j, i) * GMG(l, k) - kronecker(bl1, bl2) * GM(l, i) * MG(j, k));
      }
  }

  void M3ph_tau::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z           = mpi_all_reduce(Z, comm);
    M3ph_tau_   = mpi_all_reduce(M3ph_tau_, comm);
    double dtau = params.beta / (params.n_tau_M3 - 1);
    M3ph_tau_   = M3ph_tau_ / (Z * dtau * dtau);
  }

} // namespace triqs_ctint::measures
