#include "./f4ph_loc_iw.hpp"
#include <cmath>

namespace triqs_ctint::measures {
  f4ph_loc_iw::f4ph_loc_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, g_iw_t const &GinvG0_iw_,
                           g_iw_t const &G0Ginv_iw_)
     : params(params_), GinvG0_iw(GinvG0_iw_), G0Ginv_iw(G0Ginv_iw_), qmc_config(qmc_config_), buf_arrarr(params_.n_blocks()) {

    // Construct Matsubara mesh for temporary Matrix
    auto iw_mesh_large = mesh::imfreq{params.beta, Fermion, params.n_iW_M4 + params.n_iw_M4 + 1};
    auto M_mesh        = iw_mesh_large * iw_mesh_large;

    // Initialize intermediate scattering matrix
    M = block_gf{M_mesh, params.gf_struct};

    // Initialize buffer for Ginv * G0 * M * G0 * Ginv
    auto m_bl = gf<prod<imfreq, imfreq>, tensor_valued<1>>{M_mesh, {M[0].target_shape()[0]}};
    m         = make_block_gf(params.n_blocks(), m_bl);

    // Construct Matsubara mesh
    auto iW_mesh          = mesh::imfreq{params.beta, Boson, params.n_iW_M4};
    auto iw_mesh          = mesh::imfreq{params.beta, Fermion, params.n_iw_M4};
    auto f4ph_loc_iw_mesh = iW_mesh * iw_mesh * iw_mesh;

    // Init measurement container and capture view
    auto f4ph_loc_iw_bl  = gf<prod<imfreq, imfreq, imfreq>, tensor_valued<1>>{f4ph_loc_iw_mesh, {M[0].target_shape()[0]}};
    results->f4ph_loc_iw = make_block2_gf(params.n_blocks(), params.n_blocks(), f4ph_loc_iw_bl);
    f4ph_loc_iw_.rebind(results->f4ph_loc_iw.value());
    f4ph_loc_iw_() = 0;

    // Create nfft buffers
    for (int bl : range(params.n_blocks())) {
      auto init_target_func = [&](int i, int j) {
        return nfft_buf_t<2>{slice_target_to_scalar(M[bl], i, j).data(), params.nfft_buf_size, params.beta};
      };
      buf_arrarr(bl) = array_adapter{M[bl].target_shape(), init_target_func};
    }
  }

  void f4ph_loc_iw::accumulate(mc_weight_t sign) {
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
      auto const &L_bl = GinvG0_iw[bl];
      auto const &R_bl = G0Ginv_iw[bl];

      for (auto [iw1, iw2] : m_bl.mesh())
        for (int i : range(bl_size))
          for (int j : range(bl_size))
            for (int k : range(bl_size)) { m_bl[iw1, iw2](i) += L_bl[iw1.value()](i, j) * M_bl[iw1, iw2](j, k) * R_bl[iw2.value()](k, i); }
    }

    // Calculate f4ph_loc
    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {
        int bl_size    = m[bl1].target_shape()[0];
        auto const &m1 = m[bl1];
        auto const &m2 = m[bl2];
        auto &f4_loc   = f4ph_loc_iw_(bl1, bl2);

        for (auto [iW, iw, iwp] : f4_loc.mesh())
          for (int i : range(bl_size)) {
            f4_loc[iW, iw, iwp](i) += sign * m1[iW + iw, iw.value()](i) * m2[iwp.value(), iW + iwp](i);
            if (bl1 == bl2) { f4_loc[iW, iw, iwp](i) -= sign * m1[iwp.value(), iw.value()](i) * m2[iW + iw, iW + iwp](i); }
          }
      }
  }

  void f4ph_loc_iw::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z            = mpi::all_reduce(Z, comm);
    f4ph_loc_iw_ = mpi::all_reduce(f4ph_loc_iw_, comm);
    f4ph_loc_iw_ = f4ph_loc_iw_ / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
