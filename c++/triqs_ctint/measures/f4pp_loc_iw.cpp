#include "./f4pp_loc_iw.hpp"
#include <cmath>

namespace triqs_ctint::measures {
  f4pp_loc_iw::f4pp_loc_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_iw_t const &GinvG01_iw_,
                           g_iw_t const &GinvG02_iw_)
     : params(params_), GinvG01_iw(GinvG01_iw_), GinvG02_iw(GinvG02_iw_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh for temporary Matrix
    mesh::imfreq iw_mesh_large{params.beta, Fermion, params.n_iW_M4 + params.n_iw_M4 + 1};
    mesh::prod<imfreq, imfreq> M_mesh{iw_mesh_large, iw_mesh_large};

    // Initialize intermediate scattering matrix
    M = block_gf{M_mesh, params.gf_struct};

    // Initialize buffer for Ginv * G0 * M * G0 * Ginv
    auto m_bl = gf<prod<imfreq, imfreq>, tensor_valued<1>>{M_mesh, {M[0].target_shape()[0]}};
    m         = make_block_gf(params.n_blocks(), m_bl);

    // Construct Matsubara mesh
    mesh::imfreq iW_mesh{params.beta, Boson, params.n_iW_M4};
    mesh::imfreq iw_mesh{params.beta, Fermion, params.n_iw_M4};
    mesh::prod<imfreq, imfreq, imfreq> f4pp_loc_iw_mesh{iW_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    auto f4pp_loc_iw_bl  = gf<prod<imfreq, imfreq, imfreq>, tensor_valued<1>>{f4pp_loc_iw_mesh, {M[0].target_shape()[0]}};
    results->f4pp_loc_iw = make_block2_gf(params.n_blocks(), params.n_blocks(), f4pp_loc_iw_bl);
    f4pp_loc_iw_.rebind(results->f4pp_loc_iw.value());
    f4pp_loc_iw_() = 0;

    // Create nfft buffers
    for (int bl : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft_buf_t<2>{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta};
      };
      buf_arrarr(bl) = array_adapter{M[bl].target_shape(), init_target_func};
    }
  }

  void f4pp_loc_iw::accumulate(mc_weight_t sign) {
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

    // Calculate Ginv * G0 * M * G0 * Ginv
    m() = 0;

    for (int bl : range(params.n_blocks())) {
      int bl_size      = m[bl].target_shape()[0];
      auto &m_bl       = m[bl];
      auto const &M_bl = M[bl];
      auto const &L_bl = GinvG01_iw[bl];
      auto const &R_bl = GinvG02_iw[bl];

      for (auto iw1 : std::get<0>(m_bl.mesh()))
        for (auto iw2 : std::get<1>(m_bl.mesh()))
          for (int i : range(bl_size))
            for (int j : range(bl_size))
              for (int k : range(bl_size)) { m_bl[iw1, iw2](i) += L_bl(iw1)(i, j) * M_bl[iw1, iw2](j, k) * R_bl(iw2)(k, i); }
    }

    // Calculate f4pp_loc
    auto const &iW_mesh = std::get<0>(f4pp_loc_iw_(0, 0).mesh());
    auto const &iw_mesh = std::get<1>(f4pp_loc_iw_(0, 0).mesh());

    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {
        int bl_size    = m[bl1].target_shape()[0];
        auto const &m1 = m[bl1];
        auto const &m2 = m[bl2];
        auto &f4_loc   = f4pp_loc_iw_(bl1, bl2);

        for (auto iW : iW_mesh)
          for (auto iw : iw_mesh)
            for (auto iwp : iw_mesh)
              for (int i : range(bl_size)) {
                f4_loc[iW, iw, iwp](i) += sign * m1(iW - iwp, iw)(i) * m2(iwp, iW - iw)(i);
                if (bl1 == bl2) { f4_loc[iW, iw, iwp](i) -= sign * m1(iwp, iw)(i) * m2(iW - iwp, iW - iw)(i); }
              }
      }
  }

  void f4pp_loc_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z            = mpi::all_reduce(Z, comm);
    f4pp_loc_iw_ = mpi::all_reduce(f4pp_loc_iw_, comm);
    f4pp_loc_iw_ = f4pp_loc_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
