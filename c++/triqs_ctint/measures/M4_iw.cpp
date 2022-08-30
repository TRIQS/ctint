#include "./M4_iw.hpp"
#include <cmath>

namespace triqs_ctint::measures {

  M4_iw::M4_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    gf_mesh<imfreq> iw_mesh{params.beta, Fermion, params.n_iw_M4};
    gf_mesh<prod<imfreq, imfreq, imfreq>> M4_iw_mesh{iw_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M4_iw = make_block2_gf(M4_iw_mesh, params.gf_struct);
    M4_iw_.rebind(results->M4_iw.value());
    for (auto &&[bl1, bl2] : params.block_pairs_indices_M4()) {
      M4_iw_(bl1, bl2)() = 0;
    }

    // Construct Matsubara mesh for temporary Matrix
    gf_mesh<imfreq> iw_mesh_large{params.beta, Fermion, 3 * params.n_iw_M4};
    gf_mesh<prod<imfreq, imfreq>> M_mesh{iw_mesh_large, iw_mesh};

    // Initialize intermediate scattering matrix
    M = block_gf{M_mesh, params.gf_struct};

    // Create nfft buffers
    for (int bl : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft_buf_t<2>{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta};
      };
      buf_arrarr(bl) = array_adapter{M[bl].target_shape(), init_target_func};
    }
  }

  void M4_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    M() = 0;
    for (int bl : range(params.n_blocks()))
      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[bl],
               [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) { // Care for negative frequency in c transform (for M-objects)
                 buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({double(cdag_j.tau), params.beta - double(c_i.tau)}, -Ginv_ji);
               })
        ;
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    for (auto &&[bl1, bl2] : params.block_pairs_indices_M4()) {
      int bl1_size   = M[bl1].target_shape()[0];
      int bl2_size   = M[bl2].target_shape()[0];
      auto const &M1 = M[bl1];
      auto const &M2 = M[bl2];
      auto &M4       = M4_iw_(bl1, bl2);
      auto const &iw_mesh = std::get<0>(M4.mesh());
      for (int i : range(bl1_size))
        for (int j : range(bl1_size))
          for (int k : range(bl2_size))
            for (int l : range(bl2_size))
              for (auto const &iw1 : iw_mesh)
                for (auto const &iw2 : iw_mesh)
                  for (auto const &iw3 : iw_mesh) {
                    gf_mesh<imfreq>::mesh_point_t iw4{iw_mesh, iw1.index() + iw3.index() - iw2.index()};
                    M4[{iw1, iw2, iw3}](i, j, k, l) +=
                       sign * (M1[{iw2, iw1}](j, i) * M2[{iw4, iw3}](l, k) - kronecker(bl1, bl2) * M1[{iw4, iw1}](l, i) * M2[{iw2, iw3}](j, k));
                  }
      }
  }

  void M4_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z      = mpi::all_reduce(Z, comm);
    for (auto &&[bl1, bl2] : params.block_pairs_indices_M4()) {
      M4_iw_(bl1, bl2) = mpi::all_reduce(M4_iw_(bl1, bl2), comm);
      M4_iw_(bl1, bl2) = M4_iw_(bl1, bl2) / (Z * params.beta);
    }
  }

} // namespace triqs_ctint::measures
