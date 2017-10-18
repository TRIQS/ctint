#include "./qmc_config.hpp"

namespace triqs_ctint {

  qmc_config_t::qmc_config_t(params_t const &params, block_gf<imtime, matrix_valued> const &G0_tau) {

    // Build the determinants. Explicit value '1000' reserves matrix memory (c.f. det doc)
    for (int i = 0; i < G0_tau.size(); ++i) dets.emplace_back(G0hat_t{G0_tau[i], params.alpha[i]}, 1000);
  }

} // namespace triqs_ctint
