#include "./M_iw_dlr.hpp"
#include "triqs_ctint/types.hpp"
#include <mpi/generic_communication.hpp>
#include <nda/basic_functions.hpp>
#include <nda/concepts.hpp>
#include <triqs/gfs/block/block_gf.hpp>
#include <triqs/gfs/block/factories.hpp>
#include <triqs/mesh/utils.hpp>
#include <triqs/utility/first_include.hpp>

namespace triqs_ctint::measures {

  M_iw_dlr::M_iw_dlr(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : params(params_),
       qmc_config(qmc_config_),
       tau_samples(results->tau_samples),
       weight_samples(results->weight_samples),
       M_iw_dlr_(results->M_iw_dlr) {

    M_iw_dlr_              = block_gf{mesh::dlr_imfreq{params.beta, Fermion, params.w_max, params.eps}, params.gf_struct};
    auto n_matrix_elements = 0;
    for (auto [bl, bl_size] : params.gf_struct) n_matrix_elements += bl_size * bl_size;

    for (auto [bl, bl_size] : params.gf_struct) {
      tau_samples.emplace_back(bl_size, bl_size);
      weight_samples.emplace_back(bl_size, bl_size);
      for (int i = 0; i < bl_size; i++) {
        for (int j = 0; j < bl_size; j++) {
          tau_samples.back()(i, j).reserve(params.sample_buffer_size / n_matrix_elements);
          weight_samples.back()(i, j).reserve(params.sample_buffer_size / n_matrix_elements);
        }
      }
    }
  }

  void M_iw_dlr::accumulate(mc_weight_t sign) {
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

          if (ws_ij.size() == ws_ij.capacity()) { convert_samples_to_dlr(); }

          // Absolute time-difference tau of the index pair
          auto [s, dtau] = cyclic_difference(cdag_j.tau, c_i.tau);

          ts_ij.emplace_back(dtau);
          ws_ij.emplace_back(Ginv * s * sign);
        }
      });
    }
  }

  void M_iw_dlr::collect_results(mpi::communicator const &comm) {

    convert_samples_to_dlr();

    Z = mpi::all_reduce(Z, comm);

    // Collect results and normalize
    mpi::all_reduce_in_place(M_iw_dlr_, comm);
    M_iw_dlr_ /= (-Z * params.beta);
  }

  void M_iw_dlr::convert_samples_to_dlr() {
    auto const &dlr_mesh = M_iw_dlr_[0].mesh();
    for (auto bl : range(params.n_blocks())) {
      auto bl_size = params.gf_struct[bl].second;
      for (auto i : range(bl_size)) {
        for (auto j : range(bl_size)) {
          auto ws  = nda::array_view<dcomplex, 1>(weight_samples[bl](i, j));
          auto tmp = nda::array_view<double, 1>(tau_samples[bl](i, j));
          auto ts  = nda::array<dcomplex, 1>(tmp);
          for (auto iw : dlr_mesh) { M_iw_dlr_[bl][iw](i, j) += nda::linalg::dot(ws, nda::exp(iw * ts)); }

          // Clear the samples after processing
          weight_samples[bl](i, j).clear();
          tau_samples[bl](i, j).clear();
        }
      }
    }
  }

} // namespace triqs_ctint::measures
