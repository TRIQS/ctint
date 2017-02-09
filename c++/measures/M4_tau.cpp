#include "./M4_tau.hpp"

namespace triqs_ctint::measures {

  M4_tau::M4_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

     // Construct imaginary time mesh
    gf_mesh<imtime> time_mesh{params.beta, Fermion, params.n_tau_M4};
    gf_mesh<cartesian_product<imtime, imtime, imtime>> M4_tau_mesh{time_mesh, time_mesh, time_mesh};

    // Init measurement container and capture view
    std::vector<std::vector<block_type>> v_tau;
    for (auto const &bl1 : params.gf_struct) {
      std::vector<block_type> temp;
      for (auto const &bl2 : params.gf_struct)
        temp.emplace_back(M4_tau_mesh, make_shape(bl1.second.size(), bl1.second.size(), bl2.second.size(), bl2.second.size()));
      v_tau.emplace_back(std::move(temp));
    }

    results->M4_tau = make_block2_gf(params.block_names(), params.block_names(), v_tau);
    M4_tau_.rebind(*results->M4_tau);
    M4_tau_() = 0;
  }

  void M4_tau::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    // TODO Check: Mbar calculation notes?
    for (int b1 = 0; b1 < params.n_blocks(); ++b1)
      foreach (qmc_config.dets[b1], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv1) {
        for (int b2 = 0; b2 < params.n_blocks(); ++b2)
          foreach (qmc_config.dets[b2], [&](c_t const &c_k, cdag_t const &cdag_l, auto const &Ginv2) {
            int factor  = 1;
            double tau1 = cyclic_difference(cdag_j.tau, c_i.tau);
            if (cdag_j.tau < c_i.tau) factor *= -1;
            double tau2 = cyclic_difference(c_k.tau, cdag_j.tau); //bosonic time, no sign
            double tau3 = cyclic_difference(cdag_l.tau, c_k.tau);
            if (cdag_l.tau < c_k.tau) factor *= -1;

            auto M = M4_tau_(b1, b2)[closest_mesh_pt(tau1, tau2, tau3)]; // deduces to array_proxy<...>

            M(c_i.u, cdag_j.u, c_k.u, cdag_l.u) += Ginv1 * Ginv2 * factor * sign;
            if (b1 == b2) M(c_i.u, cdag_l.u, c_k.u, cdag_j.u) -= Ginv1 * Ginv2 * factor * sign;
          })
            ;
      })
        ;
  }

  void M4_tau::collect_results(triqs::mpi::communicator const &comm) {
    // Collect results and normalize
    Z           = mpi_all_reduce(Z, comm);
    M4_tau_     = mpi_all_reduce(M4_tau_, comm);
    double dtau = params.beta / (params.n_tau_M4 - 1);
    //M4_tau_     = M4_tau_ / (Z * dtau * dtau * dtau * params.beta); FIXME
    for (auto &M : M4_tau_) M = M / (Z * dtau * dtau * dtau * params.beta);

    // Account for edge bins beeing smaller
    auto _ = var_t{}; 
    int n = params.n_tau_M4 - 1;
    for (auto &M : M4_tau_) {
       M[0][_][_] *= 2.0; 
       M[_][0][_] *= 2.0; 
       M[_][_][0] *= 2.0; 
       M[n][_][_] *= 2.0; 
       M[_][n][_] *= 2.0; 
       M[_][_][n] *= 2.0; 
    }
  }

} // namespace triqs_ctint::measures
