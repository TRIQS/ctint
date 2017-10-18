#include "./M3pp_iw.hpp"

namespace triqs_ctint::measures {

  M3pp_iw::M3pp_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, block_gf<imtime, matrix_valued> const &G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(G0_tau_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    gf_mesh<imfreq> iw_mesh{params.beta, Fermion, params.n_iw_M3};
    gf_mesh<cartesian_product<imfreq, imfreq>> M3pp_iw_mesh{iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M3pp_iw_nfft = make_block2_gf(M3pp_iw_mesh, params.gf_struct);
    M3pp_iw_.rebind(*results->M3pp_iw_nfft);
    M3pp_iw_() = 0;

    // Initialize intermediate scattering matrix
    GM = make_block_gf(iw_mesh, params.gf_struct);

    // Create nfft buffers
    for (int b : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft_buf_t<1>{slice_target_to_scalar(GM[b], i, j).data(), params.nfft_buf_size, params.beta};
      };
      buf_arrarr(b) = array<nfft_buf_t<1>, 2>{GM[b].target_shape(), init_target_func};
    }
  }

  void M3pp_iw::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    GM() = 0;
    for (int bl : range(params.n_blocks())) {
      int bl_size = GM[bl].target_shape()[0];
      for (int b_u : range(bl_size))
        //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[bl1]) // FIXME c++17
        foreach (qmc_config.dets[bl], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) {
          auto G0_bj = G0_tau[bl][closest_mesh_pt(params.beta - double(cdag_j.tau))](b_u, cdag_j.u);
          buf_arrarr(bl)(b_u, c_i.u).push_back({params.beta - double(c_i.tau)}, G0_bj * Ginv_ji);
        })
          ;
    }
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    auto const &iw_mesh = GM[0].mesh();

    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        int bl1_size    = GM[bl1].target_shape()[0];
        int bl2_size    = GM[bl2].target_shape()[0];
        auto const &GM1 = GM[bl1];
        auto const &GM2 = GM[bl2];
        auto &M3pp_iw   = M3pp_iw_(bl1, bl2);

        for (int i : range(bl1_size))
          for (int j : range(bl1_size))
            for (int k : range(bl2_size))
              for (int l : range(bl2_size))
                for (auto const &iw1 : iw_mesh)
                  for (auto const &iw3 : iw_mesh)
                    M3pp_iw[{iw1, iw3}](i, j, k, l) +=
                       sign * (GM1[iw1](j, i) * GM2[iw3](l, k) - kronecker(bl1, bl2) * GM1[iw1](l, i) * GM2[iw3](j, k));
      }
  }

  void M3pp_iw::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi_all_reduce(Z, comm);
    M3pp_iw_ = mpi_all_reduce(M3pp_iw_, comm);
    M3pp_iw_ = M3pp_iw_ / Z;
  }

} // namespace triqs_ctint::measures
