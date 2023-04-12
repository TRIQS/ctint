#include "./chi2_tau.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  template <Chan_t Chan>
  chi2_tau<Chan>::chi2_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), tau_mesh{params_.beta, Boson, params_.n_tau_chi2} {

    // Init measurement container and capture view
    if (Chan == Chan_t::PP) {
      results->chi2pp_tau = make_block2_gf(tau_mesh, params.gf_struct);
      chi2_tau_.rebind(results->chi2pp_tau.value());
    } else if (Chan == Chan_t::PH) {
      results->chi2ph_tau = make_block2_gf(tau_mesh, params.gf_struct);
      chi2_tau_.rebind(results->chi2ph_tau.value());
    }
    chi2_tau_() = 0;
  }

  template <Chan_t Chan> void chi2_tau<Chan>::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    // Calculate chi2ph
    for (int bl1 : range(params.n_blocks()))
      for (int bl2 : range(params.n_blocks())) {

        auto &det1     = qmc_config.dets[bl1];
        auto &det2     = qmc_config.dets[bl2];
        auto &chi2_tau = chi2_tau_(bl1, bl2);
        int bl1_size   = chi2_tau.target_shape()[0];
        int bl2_size   = chi2_tau.target_shape()[2];

        for (int a : range(bl1_size))
          for (int b : range(bl1_size))
            for (int c : range(bl2_size))
              for (int d : range(bl2_size))
                for (auto tau : tau_mesh) {

                  cdag_t cdag_a, cdag_c;
                  c_t c_b, c_d;

                  auto tau_point  = make_tau_t(double(tau));
                  auto taup_point = tau_point;

                  // Time-Ordering
                  tau_point.n += 3;  // tau -> tau^{+++}
                  taup_point.n += 2; // taup -> tau^{++}

                  if (tau.idx == params.n_tau_chi2 - 1) {
                    tau_point  = tau_t::get_beta_minus();
                    taup_point = tau_t::get_beta_minus_minus();
                  }

                  if constexpr (Chan == Chan_t::PP) { // Particle-particle channel

                    cdag_a = cdag_t{tau_point, a};
                    c_b    = c_t{tau_t::get_zero_plus(), b};
                    cdag_c = cdag_t{taup_point, c};
                    c_d    = c_t{tau_t::get_zero(), d};

                  } else if constexpr (Chan == Chan_t::PH) { // Particle-hole channel

                    cdag_a = cdag_t{tau_point, a};
                    c_b    = c_t{taup_point, b};
                    cdag_c = cdag_t{tau_t::get_zero_plus(), c};
                    c_d    = c_t{tau_t::get_zero(), d};
                  }

                  if (bl1 == bl2)
                    chi2_tau[tau](a, b, c, d) += sign * det1.try_insert2(0, 1, 0, 1, c_b, c_d, cdag_a, cdag_c);
                  else
                    chi2_tau[tau](a, b, c, d) += sign * det1.try_insert(0, 0, c_b, cdag_a) * det2.try_insert(0, 0, c_d, cdag_c);

                  det1.reject_last_try();
                  det2.reject_last_try();
                }
      }
  }

  template <Chan_t Chan> void chi2_tau<Chan>::collect_results(mpi::communicator const &comm) {
    // Collect results and normalize
    Z         = mpi::all_reduce(Z, comm);
    chi2_tau_ = mpi::all_reduce(chi2_tau_, comm);
    chi2_tau_ = chi2_tau_ / Z;
  }

} // namespace triqs_ctint::measures
