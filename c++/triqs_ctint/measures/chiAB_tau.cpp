#include "./chiAB_tau.hpp"

using namespace triqs::utility;

namespace triqs_ctint::measures {

  chiAB_tau::chiAB_tau(params_t const &params_, qmc_config_t &qmc_config_, container_set *results)
     : params(params_), qmc_config(qmc_config_), tau_mesh{params_.beta, Boson, params_.n_tau_chi2} {

    if (params.chi_A_vec.size() == 0 or params.chi_B_vec.size() == 0)
      TRIQS_RUNTIME_ERROR << " Empty operator vector detected in chiAB measurement \n";

    // Function that takes a bosonic operator Op = Sum_i a_i c^+(bi, ui) c(bi, vi)
    // and returns a vector<tuple> with v[i] = (b_i, a_i, (c^+(bi, ui), c(bi, vi)))
    auto get_terms = [&params_](many_body_operator &A) {
      std::vector<op_term_t> terms;
      for (auto const &term : A) {
        auto const &m = term.monomial;
        if (m.size() != 2 or !m[0].dagger or m[1].dagger)
          TRIQS_RUNTIME_ERROR << " Monomial in bosonic operator of chiAB measurement not of the proper form c^+ c \n";
        auto [bl1, i] = get_int_indices(m[0], params_.gf_struct);
        auto [bl2, j] = get_int_indices(m[1], params_.gf_struct);
        auto bl_pair  = std::make_pair(bl1, bl2);
        auto op_pair  = std::make_pair(cdag_t{tau_t::get_zero_plus(), i}, c_t{tau_t::get_zero(), j});
        terms.emplace_back(term.coef, bl_pair, op_pair);
      }
      return terms;
    };

    for (auto A : params.chi_A_vec) A_vec.emplace_back(get_terms(A));
    for (auto B : params.chi_B_vec) B_vec.emplace_back(get_terms(B));

    // Init measurement container and capture view
    results->chiAB_tau = gf<imtime>{tau_mesh, make_shape(A_vec.size(), B_vec.size())};
    chiAB_tau_.rebind(*results->chiAB_tau);
    chiAB_tau_() = 0;
  }

  void chiAB_tau::accumulate(mc_weight_t sign) {
    // Accumulate sign
    Z += sign;

    for (auto [j, B] : enumerate(B_vec))
      for (auto &[coef_B, bl_pair_B, op_pair_B] : B) {

        auto [cdag_B, c_B]       = op_pair_B;
        auto [bl_cdag_B, bl_c_B] = bl_pair_B;

        for (auto [i, A] : enumerate(A_vec))
          for (auto &[coef_A, bl_pair_A, op_pair_A] : A) {

            auto [cdag_A, c_A]       = op_pair_A;
            auto [bl_cdag_A, bl_c_A] = bl_pair_A;

            // Determine the block order in the quartet of operators
            bool is_AABB = (bl_cdag_A == bl_c_A && bl_cdag_B == bl_c_B);
            bool is_ABAB = (bl_cdag_A == bl_c_B && bl_cdag_B == bl_c_A);
            bool is_AAAA = is_AABB && is_ABAB;

            // Define the determinants to operate on
            auto &det_1 = qmc_config.dets[bl_cdag_A];
            auto &det_2 = is_AABB ? qmc_config.dets[bl_cdag_B] : qmc_config.dets[bl_c_A];

            for (auto tau : tau_mesh) {
              auto tau_point  = make_tau_t(double(tau));
              auto taup_point = tau_point;

              // Time-Ordering
              tau_point.n += 3;  // tau -> tau^{+++}
              taup_point.n += 2; // taup -> tau^{++}

              if (tau.index() == params.n_tau_chi2 - 1) {
                tau_point  = tau_t::get_beta_minus();
                taup_point = tau_t::get_beta_minus_minus();
              }

              cdag_A.tau = tau_point;
              c_A.tau    = taup_point;

	      // TODO Extend for measurement of pairing response functions
              if (is_AAAA)
                chiAB_tau_[tau](i, j) += coef_A * coef_B * sign * det_1.try_insert2(0, 1, 0, 1, c_A, c_B, cdag_A, cdag_B);
              else if (is_AABB)
                chiAB_tau_[tau](i, j) += coef_A * coef_B * sign * det_1.try_insert(0, 0, c_A, cdag_A) * det_2.try_insert(0, 0, c_B, cdag_B);
              else if (is_ABAB) // In this case we have to swap c_A and c_B, which leads to a minus sign
                chiAB_tau_[tau](i, j) -= coef_A * coef_B * sign * det_1.try_insert(0, 0, c_B, cdag_A) * det_2.try_insert(0, 0, c_A, cdag_B);
	      // Note: All other block combinations vanish

              det_1.reject_last_try();
              det_2.reject_last_try();
            }
          }
      }
  }

  void chiAB_tau::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z          = mpi_all_reduce(Z, comm);
    chiAB_tau_ = mpi_all_reduce(chiAB_tau_, comm);
    chiAB_tau_ = chiAB_tau_ / Z;
  }

} // namespace triqs_ctint::measures
