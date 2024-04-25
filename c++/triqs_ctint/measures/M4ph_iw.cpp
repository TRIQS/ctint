#include "./M4ph_iw.hpp"
#include "./iw_accumulate.hpp"

namespace triqs_ctint::measures {

  M4ph_iw::M4ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    mesh::imfreq iW_mesh{params.beta, Boson, params.n_iW_M4};
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M4};
    mesh::prod<imfreq, imfreq, imfreq> M4ph_iw_mesh{iW_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M4ph_iw = make_block2_gf(M4ph_iw_mesh, params.gf_struct);
    M4ph_iw_.rebind(results->M4ph_iw.value());
    M4ph_iw_() = 0;

    // Construct Matsubara mesh for temporary Matrix
    mesh::imfreq iw_mesh_large{params.beta, Fermion, params.n_iW_M4 + params.n_iw_M4 + 1};
    mesh::prod<imfreq, imfreq> M_mesh{iw_mesh_large, iw_mesh_large};

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

  void M4ph_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    M() = 0;
    for (int bl : range(params.n_blocks()))
      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[bl],
               [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) { // Care for negative frequency in c transform (for M-objects)
                 buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({double(cdag_j.tau), params.beta - double(c_i.tau)}, -Ginv_ji);
               });
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    for (int bl1 : range(params.n_blocks())) { // FIXME c++17 Loops
      for (int bl2 : range(params.n_blocks())) {
        auto const bl2_size = M[bl2].target_shape()[0];
        simd::iw4ph_accumulate(sign, M, M4ph_iw_, bl1, bl2, bl2_size);
      }
    }
  }

  void M4ph_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi::all_reduce(Z, comm);
    M4ph_iw_ = mpi::all_reduce(M4ph_iw_, comm);
    M4ph_iw_ = M4ph_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
