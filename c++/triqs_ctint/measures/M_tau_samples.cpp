#include "./M_tau_samples.hpp"
#include "triqs_ctint/types.hpp"
#include <mpi/generic_communication.hpp>

namespace triqs_ctint::measures {

  M_tau_samples::M_tau_samples(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    results->M_hartree = make_block_vector<M_tau_scalar_t>(params.gf_struct);
    for (auto &m : results->M_hartree.value()) M_hartree_.push_back(m);
  }

  void M_tau_samples::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    // for (auto const & [D,M] : triqs::std::zip(config.dets, results->M_tau)) 	// C++17
    auto sample= make_block_vector<sample_t>(params.gf_struct);
    for (int b = 0; b < sample.size(); ++b) {
      // Loop over every index pair (x,y) in the determinant matrix
      // for (auto const & [x,y,Ginv] : D ) 	// C++17
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {
        // Check for the equal-time case
        if (c_i.tau == cdag_j.tau) {
          M_hartree_[b](cdag_j.u, c_i.u) += Ginv * sign;
        } else {

          // Absolut time-difference tau of the index pair
          auto [s, dtau] = cyclic_difference(cdag_j.tau, c_i.tau);

          // Project tau to closest point on the binning grid
          sample[b](cdag_j.u,c_i.u) = {dtau, Ginv * s * sign};
        }
      })
        ;
    }
    M_tau_samples_.push_back(sample);
  }

  void M_tau_samples::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z           = mpi::all_reduce(Z, comm);

    double dtau = params.beta / params.n_tau;
    auto normalize = [&](auto &s){ 
      for (auto &b : s) {
        for (auto &el : b){
          el.second = el.second / (-Z*dtau*params.beta);
        }
      }
    };

    // Post-processing by Christina and Jason would be here.
    
    for (auto &m : M_hartree_) {
      m = mpi::all_reduce(m, comm);
      m = m / (-Z * params.beta);
    }
  }

} // namespace triqs_ctint::measures
