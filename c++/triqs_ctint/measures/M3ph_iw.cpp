#include "./M3ph_iw.hpp"

namespace triqs_ctint::measures {

  M3ph_iw::M3ph_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_tau_cv_t G0_tau_)
     : params(params_),
       qmc_config(qmc_config_),
       buf_arrarr(params_.n_blocks()),
       buf_arrarr_GM(params_.n_blocks()),
       buf_arrarr_MG(params_.n_blocks()),
       G0_tau(std::move(G0_tau_)) {

    // Construct Matsubara mesh
    gf_mesh<imfreq> iw_mesh{params.beta, Fermion, params.n_iw_M3};
    gf_mesh<imfreq> iW_mesh{params.beta, Boson, params.n_iW_M3};
    gf_mesh<cartesian_product<imfreq, imfreq>> M3ph_iw_mesh{iw_mesh, iW_mesh};

    // Init measurement container and capture view
    results->M3ph_iw_nfft = make_block2_gf(M3ph_iw_mesh, params.gf_struct);
    M3ph_iw_.rebind(results->M3ph_iw_nfft.value());
    M3ph_iw_() = 0;

    // Initialize intermediate scattering matrix
    gf_mesh<imfreq> iw_mesh_large{params.beta, Fermion, params.n_iw_M3 + params.n_iW_M3};
    M  = block_gf{gf_mesh<cartesian_product<imfreq, imfreq>>{iw_mesh_large, iw_mesh}, params.gf_struct};
    GM = block_gf{iw_mesh, params.gf_struct};
    MG = block_gf{iw_mesh_large, params.gf_struct};

    auto init_target_func = [&](int bl) {
      int bl_size = GM[bl].target_shape()[0];
      return array<dcomplex, 2>(bl_size, bl_size);
    };
    GMG = array<array<dcomplex, 2>, 1>{make_shape(params.n_blocks()), init_target_func};

    // Create nfft buffers
    for (int bl : range(params.n_blocks())) {
      // M
      auto init_func_M = [&](int i, int j) { return nfft_buf_t<2>{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta}; };
      buf_arrarr(bl)   = array<nfft_buf_t<2>, 2>{M[bl].target_shape(), init_func_M};

      // GM
      auto init_func_GM = [&](int i, int j) { return nfft_buf_t<1>{slice_target_to_scalar(GM[bl], i, j).data(), params.nfft_buf_size, params.beta}; };
      buf_arrarr_GM(bl) = array<nfft_buf_t<1>, 2>{GM[bl].target_shape(), init_func_GM};

      // MG
      auto init_func_MG = [&](int i, int j) { return nfft_buf_t<1>{slice_target_to_scalar(MG[bl], i, j).data(), params.nfft_buf_size, params.beta}; };
      buf_arrarr_MG(bl) = array<nfft_buf_t<1>, 2>{MG[bl].target_shape(), init_func_MG};
    }
  }

  void M3ph_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Reset intermediate scattering matrices
    for (auto &i : GMG) { i() = 0; }
    GM() = 0;
    MG() = 0;
    M()  = 0;

    double beta = params.beta;

    // Init intermediate scattering matrices
    for (int bl : range(params.n_blocks())) {
      int bl_size = GM[bl].target_shape()[0];

      //for (auto &[c_i, cdag_j, Ginv1] : qmc_config.dets[b1]) // FIXME c++17
      foreach (qmc_config.dets[bl], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv_ji) {

        auto tau_i = double(c_i.tau);
        auto tau_j = double(cdag_j.tau);

        // Fill M, Note: Minus sign from the shift of -tau_i
        buf_arrarr(bl)(cdag_j.u, c_i.u).push_back({tau_j, beta - tau_i}, -Ginv_ji);

        //Fill GMG, GM, MG
        for (int abar_u : range(bl_size)) {
          auto G0_ia = G0_tau[bl][closest_mesh_pt(c_i.tau)](c_i.u, abar_u);
          for (int b_u : range(bl_size)) {
            // Note: Minus sign from the shift of -tau_j
            auto G0_bj = -G0_tau[bl][closest_mesh_pt(beta - tau_j)](b_u, cdag_j.u);
            GMG(bl)(abar_u, b_u) += G0_bj * Ginv_ji * G0_ia;
            // Note: Minus sign from the shift of -tau_i
            buf_arrarr_GM(bl)(b_u, c_i.u).push_back({beta - tau_i}, -G0_bj * Ginv_ji);
            buf_arrarr_MG(bl)(b_u, c_i.u).push_back({tau_j}, Ginv_ji * G0_ia);
          }
        }
      })
        ;
    }

    // Flush remaining points from all buffers
    for (auto &buf_arr : buf_arrarr)
      for (auto &buf : buf_arr) buf.flush();
    for (auto &buf_arr : buf_arrarr_GM)
      for (auto &buf : buf_arr) buf.flush();
    for (auto &buf_arr : buf_arrarr_MG)
      for (auto &buf : buf_arr) buf.flush();

    auto[iw_mesh, iW_mesh] = M3ph_iw_(0, 0).mesh();

    for (int bl1 : range(params.n_blocks())) // FIXME c++17 Loops
      for (int bl2 : range(params.n_blocks())) {

        int bl1_size     = M[bl1].target_shape()[0];
        int bl2_size     = M[bl2].target_shape()[0];
        auto const &M1   = M[bl1];
        auto const &GMG2 = GMG(bl2);
        auto const &GM1  = GM[bl1];
        auto const &MG2  = MG(bl2);
        auto &M3ph_iw    = M3ph_iw_(bl1, bl2);

        for (int i : range(bl1_size))
          for (int j : range(bl1_size))
            for (int k : range(bl2_size))
              for (int l : range(bl2_size))
                for (auto const &iw : iw_mesh)
                  for (auto const &iW : iW_mesh)
                    M3ph_iw[iw, iW](i, j, k, l) +=
                       sign * (M1[iW + iw, iw](j, i) * GMG2(l, k) - kronecker(bl1, bl2) * GM1[iw](l, i) * MG2[iW + iw](j, k));
      }
  }

  void M3ph_iw::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z        = mpi_all_reduce(Z, comm);
    M3ph_iw_ = mpi_all_reduce(M3ph_iw_, comm);
    M3ph_iw_ = M3ph_iw_ / Z;
  }

} // namespace triqs_ctint::measures
