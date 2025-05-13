#include "./M_tau_samples.hpp"
#include "triqs_ctint/types.hpp"
#include <mpi/generic_communication.hpp>

namespace triqs_ctint::measures {

  M_tau_samples::M_tau_samples(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) 
  :tau_samples(results->tau_samples), weight_samples(results->weight_samples),  params(params_), qmc_config(qmc_config_) {

    for (auto [bl, bl_size] : params.gf_struct){
      tau_samples.emplace_back(nda::array<std::vector<double>,2>(bl_size, bl_size));
      weight_samples.emplace_back(nda::array<std::vector<dcomplex>,2>(bl_size, bl_size));
    }

    results->M_hartree = make_block_vector<M_tau_scalar_t>(params.gf_struct);
    for (auto &m : results->M_hartree.value()) M_hartree_.push_back(m);
  }

  void M_tau_samples::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    for (int b = 0; b < params.gf_struct.size(); ++b) {
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
          tau_samples[b](cdag_j.u, c_i.u).push_back(dtau);
          weight_samples[b](cdag_j.u, c_i.u).push_back(Ginv * s * sign);
        }
      })
        ;
    }
  }

  void M_tau_samples::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize

    Z           = mpi::all_reduce(Z, comm);

    double dtau = params.beta / params.n_tau;
    for (auto bl : range(params.gf_struct.size())){
      auto bl_size = params.gf_struct[bl].second;
      for (auto i : range(bl_size)){
        for (auto j : range(bl_size)){
          tau_samples[bl](i,j) = mpi::gather(tau_samples[bl](i,j));
          weight_samples[bl](i,j) = mpi::gather(weight_samples[bl](i,j));
        }
      }
    }

    // bl (orbital, orbital) [cheb coefficients]
    auto curlyG = std::vector<nda::matrix<std::vector<dcomplex>>>{};

    // Normalize the weights or the results?
    // nda::vector_view(weight_samples[bl](i,j)) /= (-Z*dtau*params.beta);

    // Post-processing by Christina and Jason would be here.
    
    for (auto &m : M_hartree_) {
      m = mpi::all_reduce(m, comm);
      m = m / (-Z * params.beta);
    }
  }

} // namespace triqs_ctint::measures
