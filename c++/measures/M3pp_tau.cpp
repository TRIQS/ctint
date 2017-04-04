#include "./M3pp_tau.hpp"
#include "../itertools.hpp"

namespace triqs_ctint::measures {

  M3pp_tau::M3pp_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, block_gf<imtime, matrix_valued> const &G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(G0_tau_), tau_mesh{params_.beta, Fermion, params_.n_tau_M3} {

    // Construct Matsubara mesh
    gf_mesh<cartesian_product<imtime, imtime>> M3pp_tau_mesh{tau_mesh, tau_mesh};

    // Init measurement container and capture view
    results->M3pp_tau = make_block2_gf(M3pp_tau_mesh, params.gf_struct);
    M3pp_tau_.rebind(*results->M3pp_tau);
    M3pp_tau_() = 0;
  }

  void M3pp_tau::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    // Vectors containing the binned tau-values vec[bl][i]
    std::vector<std::vector<idx_t>> c_vec, cdag_vec;

    // The two tau meshes in one dimension
    auto const &G0_tau_mesh = G0_tau[0].mesh();
    auto const &M_tau_mesh  = tau_mesh;

    // Precompute binned tau-points
    for (auto &det : qmc_config.dets) {

      auto x_to_mesh = [&M_tau_mesh](c_t const &c_i) {
        return idx_t{gf_closest_point<imtime, int>::invoke(M_tau_mesh, closest_mesh_pt(double(c_i.tau))), c_i.u};
      };

      auto y_to_mesh = [ beta = params.beta, &G0_tau_mesh ](cdag_t const &cdag_j) {
        return idx_t{gf_closest_point<imtime, int>::invoke(G0_tau_mesh, closest_mesh_pt(beta - double(cdag_j.tau))), cdag_j.u};
      };

      c_vec.push_back(make_vector_from_range(transform(det.get_x_values(), x_to_mesh)));
      cdag_vec.push_back(make_vector_from_range(transform(det.get_y_values(), y_to_mesh)));
    }

    // The intermediate scattering matrix
    std::vector<matrix<dcomplex>> GM_vec; // GM[bl](u, i)

    // Calculate intermediate scattering matrix
    for (int bl : range(params.n_blocks())) {

      auto const &cdag = cdag_vec[bl];
      auto const &det  = qmc_config.dets[bl];
      int bl_size      = G0_tau[bl].target_shape()[0];
      int det_size     = det.size();
      auto G           = matrix<dcomplex>(bl_size, det_size);

      for (int b_u : range(bl_size))
        for (int j : range(det_size)) G(b_u, j) = G0_tau[bl][cdag[j].tau_idx](b_u, cdag[j].u);

      GM_vec.push_back(G * det.inverse_matrix_values());
    }

    // Calculate M3pp
    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        auto const &GM1 = GM_vec[bl1];
        auto const &GM2 = GM_vec[bl2];
        auto const &c1  = c_vec[bl1];
        auto const &c2  = c_vec[bl2];
        int bl1_size    = G0_tau[bl1].target_shape()[0];
        int bl2_size    = G0_tau[bl2].target_shape()[0];
        int det1_size   = qmc_config.dets[bl1].size();
        int det2_size   = qmc_config.dets[bl2].size();
        auto &M3pp_tau  = M3pp_tau_(bl1, bl2);

        for (int j : range(bl1_size))
          for (int l : range(bl2_size))
            for (int i : range(det1_size))
              for (int k : range(det2_size))
                M3pp_tau[{c1[i].tau_idx, c2[k].tau_idx}](c1[i].u, j, c2[k].u, l) +=
                   sign * (GM1(j, i) * GM2(l, k) - kronecker(bl1, bl2) * GM1(l, i) * GM2(j, k));
      }
  }

  void M3pp_tau::collect_results(triqs::mpi::communicator const &comm) {

    // Collect results and normalize
    Z           = mpi_all_reduce(Z, comm);
    M3pp_tau_   = mpi_all_reduce(M3pp_tau_, comm);
    double dtau = params.beta / (params.n_tau_M3 - 1);
    M3pp_tau_   = M3pp_tau_ / (Z * dtau * dtau);
  }

} // namespace triqs_ctint::measures
