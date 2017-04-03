#include "./M_iw.hpp"

namespace triqs_ctint::measures {

  M_iw::M_iw(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    // Init Measurement Container and Capture View
    results->M_iw_nfft = make_block_gf(gf_mesh<imfreq>{params.beta, Fermion, params.n_iw}, params.gf_struct);
    M_iw_.rebind(*results->M_iw_nfft);
    M_iw_() = 0;

    // Create nfft buffers
    for (auto &b : M_iw_) {
      // Helper function to initialize array<nfft_buf_t<1>, 2>
      auto init_func = [&](int i, int j) { return triqs::utility::nfft_buf_t<1>{b.data()(range(), i, j), params.nfft_buf_size, params.beta, true}; };
      // Initialize vector of array<nfft_buf_t<1>, 2>
      buf_vec.emplace_back(b.target_shape(), init_func);
    }
  }

  void M_iw::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    // for (auto const & [D,M] : triqs::std::zip(qmc_config.dets, results->M_tau)) 	// C++17
    for (int b = 0; b < M_iw_.size(); b++) {
      // Loop over every index pair (x,y) in the determinant matrix
      // for (auto const & [x,y,Ginv] : D ) 	// C++17
      foreach (qmc_config.dets[b], [&](arg_t const &x, arg_t const &y, auto const &Ginv) {

        // Absolut time-difference tau of the index pair
        double tau = cyclic_difference(y.tau, x.tau);

        // Care for sign-change in case of tau-shift
        int factor = (x.tau > y.tau) ? -sign : sign;

        // Push {tau, f(tau)} pair into nfft buffer
        auto &buf = buf_vec[b](y.a, x.a);
        buf.push_back({tau}, Ginv * factor);
      })
        ;
    }
  }

  void M_iw::collect_results(triqs::mpi::communicator const &comm) {
    // Flush remaining points in nfft buffers
    for (auto &buf_arr : buf_vec)
      for (auto &buf : buf_arr) buf.flush();

    // Collect results and normalize
    Z     = mpi_all_reduce(Z, comm);
    M_iw_ = mpi_all_reduce(M_iw_, comm);
    M_iw_ = M_iw_ / (-Z * params.beta);
  }

} // namespace triqs_ctint::measures
