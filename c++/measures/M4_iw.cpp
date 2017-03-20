#include "./M4_iw.hpp"

namespace triqs_ctint::measures {

  M4_iw::M4_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh
    gf_mesh<imfreq> iw_mesh{params.beta, Fermion, params.n_iw_M4};
    gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> M4_iw_mesh{iw_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    results->M4_iw = make_block2_gf(M4_iw_mesh, params.gf_struct);
    M4_iw_.rebind(*results->M4_iw);
    M4_iw_() = 0;

    // Construct Matsubara mesh for temporary Matrix
    gf_mesh<imfreq> iw_mesh_large{params.beta, Fermion, 3 * params.n_iw_M4};
    gf_mesh<cartesian_product<imfreq, imfreq>> M_mesh{iw_mesh_large, iw_mesh};

    // Initialize intermediate scattering matrix
    M = make_block_gf(M_mesh, params.gf_struct);

    // Create nfft buffers
    for (int b : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft_buf_t<2>{slice_target_to_scalar(M[b], i, j).data(), params.nfft_buf_size, params.beta};
      };
      buf_arrarr(b) = array<nfft_buf_t<2>, 2>{M[b].target_shape(), init_target_func};
    }
  }

  void M4_iw::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    // Calculate intermediate scattering matrix
    M() = 0;
    for (int b : range(params.n_blocks()))
      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) { // Care for negative frequency in Cdag transform
        buf_arrarr(b)(cdag_j.u, c_i.u).push_back({params.beta - double(cdag_j.tau), double(c_i.tau)}, Ginv);
      })
        ;
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush(); // Flush remaining points from all buffers

    for (int bl1 : range(params.n_blocks())) // FIXME Block indeces for blocks working?
      for (int bl2 : range(params.n_blocks()))
        M4_iw_(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_) << M4_iw_(bl1, bl2)(iw1_, iw2_, iw3_)(i_, j_, k_, l_)
              + sign * (M(bl1)(iw2_, iw1_)(j_, i_) * M(bl2)(iw1_ + iw3_ - iw2_, iw3_)(l_, k_)
                        - kronecker(bl1, bl2) * M(bl1)(iw1_ + iw3_ - iw2_, iw1_)(l_, i_) * M(bl2)(iw2_, iw3_)(j_, k_));
  }

  void M4_iw::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z      = mpi_all_reduce(Z, comm);
    M4_iw_ = mpi_all_reduce(M4_iw_, comm);
    M4_iw_ = M4_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
