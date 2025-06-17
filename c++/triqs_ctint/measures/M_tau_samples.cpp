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

  void M_tau_samples::collect_results(mpi::communicator const &comm) {

    auto p                    = params.n_cheb_coeffs;
    int N                     = 0; // Number of samples in each block, orbital index pair
    auto map_to_cheb_interval = [a = 0, b = params.beta](auto const &x) { return nda::array<double, 1>{(2 * x - (a + b)) / (b - a)}; };
    for (auto bl : range(params.gf_struct.size())) {
      auto bl_size = params.gf_struct[bl].second;
      for (auto i : range(bl_size)) {
        for (auto j : range(bl_size)) {
          N                = tau_samples[bl](i, j).size();
          auto T_next      = nda::array<double, 1>(N);
          auto weight      = nda::array_view<dcomplex, 1>(weight_samples[bl](i, j));
          auto tau         = nda::array_view<double, 1>(tau_samples[bl](i, j));
          auto tau_prime   = map_to_cheb_interval(tau);
          auto cheb_weight = 1.0 / nda::sqrt(1.0 - tau_prime * tau_prime);
          weight *= cheb_weight;
          // Begin Chebyshev recurrence
          auto T_prev = nda::ones<double>(N);
          curlyG[bl](i, j).push_back(nda::dot(weight, T_prev) / (M_PI));
          auto T_curr = tau_prime; //0, params.beta);
          curlyG[bl](i, j).push_back(nda::dot(weight, T_curr) / (M_PI / 2));
          for (int k = 2; k < p + 1; k++) {
            T_next = 2 * map_to_cheb_interval(tau) * T_curr - T_prev;
            curlyG[bl](i, j).push_back(nda::dot(weight, T_next) / (M_PI / 2));
            T_prev = T_curr;
            T_curr = T_next;
          }
        }
      }
    }

    Z = mpi::all_reduce(Z, comm);

    // Collect results and normalize
    for (auto b : range(params.gf_struct.size())) {
      auto bl_size = params.gf_struct[b].second;
      for (auto i : range(bl_size)) {
        for (auto j : range(bl_size)) {
          tau_samples[b](i, j)    = mpi::gather(tau_samples[b](i, j));
          weight_samples[b](i, j) = mpi::gather(weight_samples[b](i, j));
          curlyG[b](i, j)         = mpi::reduce(curlyG[b](i, j), comm);
          nda::vector_view<dcomplex>{curlyG[b](i, j)} /= (-Z * params.beta);
          // curlyG[b](i, j)         = tmp/(-Z * params.beta);
        }
      }
    }

    for (auto &m : M_hartree_) {
      m = mpi::all_reduce(m, comm);
      m = m / (-Z * params.beta);
    }
  }

} // namespace triqs_ctint::measures
