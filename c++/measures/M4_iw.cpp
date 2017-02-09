#include "./M4_iw.hpp"

namespace triqs_ctint::measures {

  M4_iw::M4_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    // Construct Matsubara mesh
    gf_mesh<imfreq> iw_mesh{params.beta, Fermion, params.n_iw_M4};
    gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> M4_iw_mesh{iw_mesh, iw_mesh, iw_mesh};

    // Init measurement container and capture view
    std::vector<std::vector<block_type>> v_iw;
    for (auto const &bl1 : params.gf_struct) {
      std::vector<block_type> temp_vec;
      for (auto const &bl2 : params.gf_struct)
        temp_vec.emplace_back(M4_iw_mesh, make_shape(bl1.second.size(), bl1.second.size(), bl2.second.size(), bl2.second.size()));
      v_iw.emplace_back(std::move(temp_vec));
    }
    results->M4_iw_nfft = make_block2_gf(params.block_names(), params.block_names(), v_iw);
    M4_iw_.rebind(*results->M4_iw_nfft);
    M4_iw_() = 0;

    // Create nfft buffers
    for (int b1 = 0; b1 < params.n_blocks(); ++b1) {
      std::vector<array<triqs::utility::nfft_buf_t<3>, 4>> temp_vec;
      for (int b2 = 0; b2 < params.n_blocks(); ++b2) {
        auto init_target_func = [&](int i, int j, int k, int l) {
          return triqs::utility::nfft_buf_t<3>{slice_target_to_scalar(M4_iw_(b1, b2), i, j, k, l).data(), params.nfft_buf_size, params.beta};
        };
        temp_vec.emplace_back(M4_iw_(b1, b2).target_shape(), init_target_func);
      }
      buf_vecvec.emplace_back(std::move(temp_vec));
    }
  }

  void M4_iw::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    for (int b1 = 0; b1 < params.n_blocks(); ++b1)
      foreach (qmc_config.dets[b1], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv1) {
        for (int b2 = 0; b2 < params.n_blocks(); ++b2)
          foreach (qmc_config.dets[b2], [&](c_t const &c_k, cdag_t const &cdag_l, auto const &Ginv2) {
            int factor  = 1;
            double tau1 = cyclic_difference(cdag_j.tau, c_i.tau);
            if (cdag_j.tau < c_i.tau) factor *= -1;
            double tau2 = cyclic_difference(c_k.tau, cdag_j.tau); //bosonic time, no sign
            double tau3 = cyclic_difference(cdag_l.tau, c_k.tau);
            if (cdag_l.tau < c_k.tau) factor *= -1;
            
	    buf_vecvec[b1][b2](c_i.u, cdag_j.u, c_k.u, cdag_l.u)
	       .push_back({tau1, tau2, tau3}, Ginv1 * Ginv2 * factor * sign);
	    if (b1 == b2)
	      buf_vecvec[b1][b2](c_i.u, cdag_l.u, c_k.u, cdag_j.u)
		 .push_back({tau1, tau2, tau3}, -Ginv1 * Ginv2 * factor * sign); 
          })
            ;
      })
        ;
  }

  void M4_iw::collect_results(triqs::mpi::communicator const &comm) {
    // Flush remaining points in nfft buffers
    for (auto &buf_vec : buf_vecvec)
      for (auto &buf_arr : buf_vec)
        for (auto &buf : buf_arr) buf.flush();

    // Collect results and normalize
    Z           = mpi_all_reduce(Z, comm);
    M4_iw_      = mpi_all_reduce(M4_iw_, comm);
    double dtau = params.beta / (params.n_tau_M4 - 1);
    //M4_iw_     = M4_iw_ / (Z * params.beta); FIXME in Block2_gf
    for (auto &M : M4_iw_) M = M / (Z * params.beta);
  }

} // namespace triqs_ctint::measures
