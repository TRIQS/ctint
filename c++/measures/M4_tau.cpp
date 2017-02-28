#include "./M4_tau.hpp"

namespace triqs_ctint::measures {

  M4_tau::M4_tau(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results) : params(params_), qmc_config(qmc_config_) {

    // Construct imaginary time mesh
    gf_mesh<imtime> time_mesh{params.beta, Fermion, params.n_tau_M4};
    gf_mesh<cartesian_product<imtime, imtime, imtime>> M4_tau_mesh{time_mesh, time_mesh, time_mesh};

    // Init measurement container and capture view
    results->M4_tau = make_block2_gf(M4_tau_mesh, params.gf_struct);
    M4_tau_.rebind(*results->M4_tau);
    M4_tau_() = 0;
  }

  void M4_tau::accumulate(double sign) {
    // Accumulate sign
    Z += sign;

    for (int b1 : range(params.n_blocks()))
      foreach (qmc_config.dets[b1], [&](c_t const &c_i, cdag_t const &cdag_j, auto const &Ginv1) {
        for (int b2 : range(params.n_blocks()))
          foreach (qmc_config.dets[b2], [&](c_t const &c_k, cdag_t const &cdag_l, auto const &Ginv2) {

            auto add_to_buf = [&](auto &c_1, auto &cdag_2, auto &c_3, auto &cdag_4, double factor) {
              double tau1    = cyclic_difference(c_1.tau, cdag_4.tau);
              double tau2    = cyclic_difference(cdag_2.tau, cdag_4.tau);
              double tau3    = cyclic_difference(c_3.tau, cdag_4.tau);
              int sign_flips = int(c_1.tau < cdag_4.tau) + int(cdag_2.tau < cdag_4.tau) + int(c_3.tau < cdag_4.tau);

              //// Old Code
              //double tau1 = cyclic_difference(c_1.tau, cdag_2.tau);
              //double tau2 = cyclic_difference(cdag_2.tau, c_3.tau); // bosonic time
              //double tau3 = cyclic_difference(c_3.tau, cdag_4.tau);
              //int factor  = int(c_1.tau < cdag_2.tau) + int(c_3.tau < cdag_4.tau);

              M4_tau_(b1, b2)[closest_mesh_pt(tau1, tau2, tau3)] += Ginv1 * Ginv2 * (sign_flips % 2 ? factor : -factor);
            };

            add_to_buf(c_i, cdag_j, c_k, cdag_l, sign);
            if (b1 == b2) add_to_buf(c_i, cdag_l, c_k, cdag_j, -sign);
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
    M4_tau_     = M4_tau_ / (Z * dtau * dtau * dtau * params.beta);

    // Account for edge bins beeing smaller
    auto _ = var_t{};
    int n  = params.n_tau_M4 - 1;
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
