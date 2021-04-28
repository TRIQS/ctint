#include "./G2_fluct_diag_tau.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  G2_fluct_diag_tau::G2_fluct_diag_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), tau_mesh{params_.beta, Fermion, params_.n_tau_G2_fluct_diag} {

    // Init measurement container and capture view
    results->G2_fluct_diag_tau = make_block2_gf(tau_mesh, params.gf_struct);
    G2_fluct_diag_tau_.rebind(results->G2_fluct_diag_tau.value());

    G2_fluct_diag_tau_() = 0;
  }

  void G2_fluct_diag_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate chi2ph
    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        auto &det1     = qmc_config.dets[bl1];
        auto &det2     = qmc_config.dets[bl2];
        auto &G2_fluct_diag_tau = G2_fluct_diag_tau_(bl1, bl2);
        int bl1_size   = G2_fluct_diag_tau.target_shape()[0];
        int bl2_size   = G2_fluct_diag_tau.target_shape()[2];

        for (int a : range(bl1_size))
          for (int b : range(bl1_size))
            for (int c : range(bl2_size))
              for (auto tau : tau_mesh) {

                cdag_t cdag_a, cdag_c;
                c_t c_b, c_d;

                auto tau_point  = make_tau_t(double(tau));

                if (tau.index() == params.n_tau_G2_fluct_diag - 1) {
                  tau_point  = tau_t::get_beta_minus();
                }

                cdag_a = cdag_t{tau_point, a};
                c_b    = c_t{tau_t::get_zero_plus_plus(), b};
                cdag_c = cdag_t{tau_t::get_zero_plus(), c};
                c_d    = c_t{tau_t::get_zero(), c};

                if (bl1 == bl2)
                  G2_fluct_diag_tau[tau](a, b, c, c) += sign * det1.try_insert2(0, 1, 0, 1, c_b, c_d, cdag_a, cdag_c);
                else
                  G2_fluct_diag_tau[tau](a, b, c, c) += sign * det1.try_insert(0, 0, c_b, cdag_a) * det2.try_insert(0, 0, c_d, cdag_c);

                det1.reject_last_try();
                det2.reject_last_try();
              }
      }
  }

  void G2_fluct_diag_tau::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z         = mpi::all_reduce(Z, comm);
    G2_fluct_diag_tau_ = mpi::all_reduce(G2_fluct_diag_tau_, comm);
    G2_fluct_diag_tau_ = G2_fluct_diag_tau_ / Z;
  }

} // namespace triqs_ctint::measures
