#include "./F_tau.hpp"

namespace triqs_ctint::measures {

  F_tau::F_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results, block_gf<imtime, matrix_valued> const &G0_tau_)
     : params(params_), qmc_config(qmc_config_), G0_tau(G0_tau_) {

    // Init Measurement Container and Capture View
    results->F_tau = make_block_gf(gf_mesh<imtime>{params.beta, Fermion, params.n_tau}, params.gf_struct);
    F_tau_.rebind(*results->F_tau);
    F_tau_() = 0;
  }

  void F_tau::accumulate(double sign) {
    // sign accumulation
    Z += sign;

    // Loop over blocks
    // for (auto const & [D,M] : triqs::std::zip(config.dets, results->M_tau)) 	// C++17
    for (int b = 0; b < F_tau_.size(); ++b) {
      // Loop over every index pair (x,y) in the determinant matrix
      // for (auto const & [x,y,Ginv] : D ) 	// C++17
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {

        for (int v = 0; v < F_tau_[b].target_shape()[0]; v++)
          F_tau_[b][closest_mesh_pt(cdag_j.tau)](cdag_j.u, v) += Ginv * sign * G0_tau[b][closest_mesh_pt(c_i.tau)](c_i.u, v);
      })
        ;
    }
  }

  // finalize the Green's function construction
  void F_tau::collect_results(triqs::mpi::communicator const &comm) {

    Z           = mpi_all_reduce(Z, comm);
    F_tau_      = mpi_all_reduce(F_tau_, comm);
    double dtau = F_tau_[0].mesh().delta();
    F_tau_      = F_tau_ / (-Z * dtau);

    // Multiply first and last bin by two
    for (auto &F : F_tau_) {
      F[0] *= 2.0;
      F[params.n_tau - 1] *= 2.0;
    }
  }

} // namespace triqs_ctint::measures
