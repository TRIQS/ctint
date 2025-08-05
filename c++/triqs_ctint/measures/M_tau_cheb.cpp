#include "./M_tau_cheb.hpp"
#include "triqs_ctint/types.hpp"
#include <mpi/generic_communication.hpp>
#include <nda/basic_functions.hpp>
#include <nda/concepts.hpp>
#include <triqs/gfs/block/block_gf.hpp>
#include <triqs/gfs/block/factories.hpp>
#include <triqs/mesh/utils.hpp>
#include <triqs/stat/accumulator.hpp>
#include <triqs/utility/first_include.hpp>
#include <format>

namespace triqs_ctint::measures {

  M_tau_cheb::M_tau_cheb(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), tau_samples(results->tau_samples), weight_samples(results->weight_samples), curlyG(results->curlyG) {

    auto n_matrix_elements = 0;
    for (auto [bl, bl_size] : params.gf_struct) n_matrix_elements += bl_size * bl_size;
    //fmt::print("buffer size per matrix element: {}\n", params.sample_buffer_size / n_matrix_elements);
    for (auto [bl, bl_size] : params.gf_struct) {
      tau_samples.emplace_back(bl_size, bl_size);
      weight_samples.emplace_back(bl_size, bl_size);
      curlyG.emplace_back(nda::zeros<dcomplex>(bl_size, bl_size, params.n_cheb_coeffs));
      auto_corr_times.emplace_back(nda::zeros<double>(bl_size, bl_size, params.n_cheb_coeffs));
      curlyG_acc.emplace_back(nda::array<triqs::stat::accumulator<dcomplex>, 3>(bl_size, bl_size, params.n_cheb_coeffs));
      curlyG_acc.back() = triqs::stat::accumulator<dcomplex>{dcomplex(0.0), -1};
      for (int i = 0; i < bl_size; i++) {
        for (int j = 0; j < bl_size; j++) {
          tau_samples.back()(i, j).reserve(params.sample_buffer_size / n_matrix_elements);
          weight_samples.back()(i, j).reserve(params.sample_buffer_size / n_matrix_elements);
        }
      }
    }
  }

  void M_tau_cheb::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Loop over blocks
    for (int bl = 0; bl < params.n_blocks(); ++bl) {

      // Loop over every index pair (x,y) in the determinant matrix[b]
      foreach (qmc_config.dets[bl], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv) {
        // Check for the equal-time case
        if (c_i.tau != cdag_j.tau) { // Ignore M_hartree Contributions
          auto &ts_ij = tau_samples[bl](cdag_j.u, c_i.u);
          auto &ws_ij = weight_samples[bl](cdag_j.u, c_i.u);

          if (ws_ij.size() == ws_ij.capacity()) { convert_samples_to_coeffs(); }

          // Absolute time-difference tau of the index pair
          auto [s, dtau] = cyclic_difference(cdag_j.tau, c_i.tau);

          // Store the samples
          ts_ij.emplace_back(dtau);
          ws_ij.emplace_back(Ginv * s * sign);
        }
      });
    }
  }

  void M_tau_cheb::collect_results(mpi::communicator const &comm) {

    convert_samples_to_coeffs();

    Z = mpi::all_reduce(Z, comm);

    // Collect results and normalize
    mpi::all_reduce_in_place(curlyG, comm);
    for (auto b : range(params.n_blocks())) curlyG[b] /= (-Z * params.beta);

    for (auto bl : range(params.n_blocks())) {
      auto bl_size = params.gf_struct[bl].second;
      for (auto i : range(bl_size)) {
        for (auto j : range(bl_size)) {
          for (auto p : range(params.n_cheb_coeffs)) {
            auto &acc           = curlyG_acc[bl](i, j, p);
            auto [errs, counts] = acc.log_bin_errors_all_reduce(comm);
            // Compute the auto-correlation time for each coefficient
            auto_corr_times[bl](i, j, p) = triqs::stat::tau_estimate_from_errors(errs[int(0.7 * errs.size())], errs[0]);
            // Reset the accumulator
            acc = triqs::stat::accumulator<dcomplex>{dcomplex(0.0), -1};
            fmt::print("Block {}, Coeff {}: Auto-correlation time = {:.4f}\n", bl, p, auto_corr_times[bl](i, j, p));
          }
        }
      }
    }
  }

  void M_tau_cheb::convert_samples_to_coeffs() {

    auto p                    = params.n_cheb_coeffs;
    auto map_to_cheb_interval = [a = 0, b = params.beta](auto const &x) { return nda::array<double, 1>{(2 * x - (a + b)) / (b - a)}; };
    for (auto bl : range(params.n_blocks())) {
      auto bl_size = params.gf_struct[bl].second;
      for (auto i : range(bl_size)) {
        for (auto j : range(bl_size)) {
          int N            = tau_samples[bl](i, j).size();
          auto weight      = nda::array_view<dcomplex, 1>(weight_samples[bl](i, j));
          auto tau         = nda::array_view<double, 1>(tau_samples[bl](i, j));
          auto tau_prime   = map_to_cheb_interval(tau);
          auto cheb_weight = 1.0 / nda::sqrt(1.0 - tau_prime * tau_prime);
          weight *= cheb_weight / M_PI;

          // Prepare real and imag weight arrays
          auto w_real = make_regular(real(weight));
          auto w_imag = make_regular(imag(weight));

          // Begin Chebyshev recurrence
          auto T_prev = nda::ones<double>(N);
          auto T_curr = tau_prime;

          auto dot_prod_w = [&](auto const &v) { return nda::blas::dot(w_real, v) + 1i * nda::blas::dot(w_imag, v); };
          curlyG[bl](i, j, 0) += dot_prod_w(T_prev);
          curlyG_acc[bl](i, j, 0) << dot_prod_w(T_prev);
          curlyG[bl](i, j, 1) += (2 * dot_prod_w(T_curr));
          curlyG_acc[bl](i, j, 1) << (2 * dot_prod_w(T_curr));

          auto two_tau_prime = make_regular(2.0 * tau_prime);
          for (int k = 2; k < p; k++) {
            T_prev = two_tau_prime * T_curr - T_prev;
            curlyG[bl](i, j, k) += 2 * dot_prod_w(T_prev);
            curlyG_acc[bl](i, j, k) << (2 * dot_prod_w(T_prev));
            std::swap(T_prev, T_curr); // Swap T_prev and T_curr
          }

          // Clear the samples after processing
          weight_samples[bl](i, j).clear();
          tau_samples[bl](i, j).clear();
        }
      }
    }
  }

} // namespace triqs_ctint::measures
