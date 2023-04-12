#include "./M_iw.hpp"

namespace triqs_ctint::measures {

  M_iw::M_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    // Init measurement container and capture view
    results->M_iw_nfft = block_gf<imfreq>{{params.beta, Fermion, params.n_iw}, params.gf_struct};
    M_iw_.rebind(results->M_iw_nfft.value());
    M_iw_() = 0;

    // Create nfft buffers
    for (auto &b : M_iw_) {
      // Helper function to initialize array<nfft_buf_t<1>, 2>
      auto init_func = [&](int i, int j) { return nfft_buf_t<1>{b.data()(range::all, i, j), params.nfft_buf_size, params.beta}; };
      // Initialize vector of array<nfft_buf_t<1>, 2>
      buf_vec.emplace_back(array_adapter{b.target_shape(), init_func});
    }
  }

  void M_iw::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    // for (auto const & [D,M] : triqs::std::zip(qmc_config.dets, results->M_tau)) 	// C++17
    for (int b = 0; b < M_iw_.size(); ++b) {
      // Loop over every index pair (x,y) in the determinant matrix
      // for (auto const & [x,y,Ginv] : D ) 	// C++17
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {

        // Absolut time-difference tau of the index pair
        auto [s, dtau] = cyclic_difference(cdag_j.tau, c_i.tau);

        // Push {tau, f(tau)} pair into nfft buffer
        auto &buf = buf_vec[b](cdag_j.u, c_i.u);
        buf.push_back({dtau}, Ginv * s * sign);
      })
        ;
    }
  }

  void M_iw::collect_results(mpi::communicator const &comm) {
    // Flush remaining points in nfft buffers
    for (auto &buf_arr : buf_vec)
      for (auto &buf : buf_arr) buf.flush();

    // Collect results and normalize
    Z     = mpi::all_reduce(Z, comm);
    M_iw_ = mpi::all_reduce(M_iw_, comm);
    M_iw_ = M_iw_ / (-Z * params.beta);
  }

} // namespace triqs_ctint::measures
