#include "./M_tau_samples.hpp"
#include "triqs_ctint/types.hpp"
#include <mpi/generic_communication.hpp>
#include <nda/basic_functions.hpp>
#include <nda/concepts.hpp>

namespace triqs_ctint::measures {

  M_tau_samples::M_tau_samples(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), tau_samples(results->tau_samples), weight_samples(results->weight_samples), curlyG(results->curlyG) {

    for (auto [bl, bl_size] : params.gf_struct) {
      tau_samples.emplace_back(bl_size, bl_size);
      weight_samples.emplace_back(bl_size, bl_size);
      curlyG.emplace_back(nda::array<std::vector<dcomplex>, 2>(bl_size, bl_size));
    }

    results->M_hartree = make_block_vector<M_tau_scalar_t>(params.gf_struct);
    for (auto &m : results->M_hartree.value()) M_hartree_.push_back(m);
  }

  void M_tau_samples::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    for (int b = 0; b < params.gf_struct.size(); ++b) {

      // Loop over every index pair (x,y) in the determinant matrix[b]
      foreach (qmc_config.dets[b], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {
        // Check for the equal-time case
        if (c_i.tau == cdag_j.tau) {
          M_hartree_[b](cdag_j.u, c_i.u) += Ginv * sign;
        } else {
          // Absolute time-difference tau of the index pair
          auto [s, dtau] = cyclic_difference(cdag_j.tau, c_i.tau);
          tau_samples[b](cdag_j.u, c_i.u).emplace_back(dtau);
          weight_samples[b](cdag_j.u, c_i.u).emplace_back(Ginv * s * sign);
        }
      });
    }
  }

  nda::array<double, 1> map_to_cheb_interval(nda::array_const_view<double, 1> x, double a, double b) { return (2 * x - (a + b)) / (b - a); }
  void M_tau_samples::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize

    Z = mpi::all_reduce(Z, comm);

    double dtau = params.beta / params.n_tau;
    for (auto b : range(params.gf_struct.size())) {
      auto bl_size = params.gf_struct[b].second;
      for (auto i : range(bl_size)) {
        for (auto j : range(bl_size)) {
          tau_samples[b](i, j)    = mpi::gather(tau_samples[b](i, j));
          weight_samples[b](i, j) = mpi::gather(weight_samples[b](i, j));
        }
      }
    }

    // bl (orbital, orbital) [cheb coefficients]
    // auto curlyG = std::vector<nda::matrix<std::vector<dcomplex>>>{};
    //attempt at initializing dimensions of curlyG:
    // for (auto [bl, bl_size] : params.gf_struct) { curlyG.emplace_back(nda::array<std::vector<dcomplex>, 2>(bl_size, bl_size)); }
    // Normalize the weights or the results?
    // nda::vector_view(weight_samples[bl](i,j)) /= (-Z*dtau*params.beta);

    // Post-processing by Christine and Jason would be here.
    // for now, hard code degree of Cheb
    int p = 100;

    for (auto bl : range(params.gf_struct.size())) {
      auto bl_size = params.gf_struct[bl].second;
      for (auto i : range(bl_size)) {
        for (auto j : range(bl_size)) {
          auto N      = tau_samples[bl](i, j).size();
          auto T_prev = nda::ones<double>(N);
          //cheb_coeffs[0] = nda::dot(weight_samples[bl](i,j),T_prev)/M_PI;
          curlyG[bl](i, j).push_back(nda::dot(nda::basic_array_view{weight_samples[bl](i, j)}, T_prev) / (M_PI));
          auto T_curr = map_to_cheb_interval(nda::basic_array_view{tau_samples[bl](i, j)}, 0, params.beta);
          curlyG[bl](i, j).push_back(nda::dot(nda::basic_array_view{weight_samples[bl](i, j)}, T_curr) / (M_PI * 2));
          //cheb_coeffs[1] = nda::dot(weight_samples[bl](i,j),T_curr)/(M_PI*2);

          for (int k = 2; k < p + 1; k++) {
            auto T_next = 2 * map_to_cheb_interval(nda::basic_array_view{tau_samples[bl](i, j)}, 0, params.beta) * T_curr - T_prev;
            curlyG[bl](i, j).push_back(nda::dot(nda::basic_array_view{weight_samples[bl](i, j)}, T_next) / (M_PI * 2));
            T_prev = T_curr;
            T_curr = T_next;
          }

          //       curlyG[bl](i,j) = cheb_coeffs;
        }
      }
    }
    for (auto &m : M_hartree_) {
      m = mpi::all_reduce(m, comm);
      m = m / (-Z * params.beta);
    }
  }

} // namespace triqs_ctint::measures
