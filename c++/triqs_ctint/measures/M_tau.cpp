#include "./M_tau.hpp"

namespace triqs_ctint::measures {

  M_tau::M_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    // Init measurement container and capture view
    results->M_tau = block_gf<imtime, M_tau_target_t>{{params.beta, Fermion, params.n_tau}, params.gf_struct};
    M_tau_.rebind(*results->M_tau);
    M_tau_() = 0;
  }

  void M_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    // for (auto const & [D,M] : triqs::std::zip(config.dets, results->M_tau)) 	// C++17
    for (int b = 0; b < M_tau_.size(); ++b) {
      // Loop over every index pair (x,y) in the determinant matrix
      // for (auto const & [x,y,Ginv] : D ) 	// C++17
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {

        // Absolut time-difference tau of the index pair
        double tau = cyclic_difference(cdag_j.tau, c_i.tau);

        // Care for sign-change in case of tau-shift
        mc_weight_t factor = (cdag_j.tau < c_i.tau) ? -sign : sign;

        // Project tau to closest point on the binning grid
        M_tau_[b][closest_mesh_pt(tau)](cdag_j.u, c_i.u) += Ginv * factor;
      })
        ;
    }
  }

  void M_tau::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z           = mpi_all_reduce(Z, comm);
    M_tau_      = mpi_all_reduce(M_tau_, comm);
    double dtau = M_tau_[0].mesh().delta();
    M_tau_      = M_tau_ / (-Z * dtau * params.beta);

    // Multiply first and last bin by two
    for (auto &M : M_tau_) {
      M[0] *= 2.0;
      M[params.n_tau - 1] *= 2.0;
    }
  }

} // namespace triqs_ctint::measures
