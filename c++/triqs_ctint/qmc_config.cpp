#include "./qmc_config.hpp"

namespace triqs_ctint {

  qmc_config_t::qmc_config_t(params_t const &params, g_tau_cv_t G0_tau) {

    // Build the determinants. Explicit value '1000' reserves matrix memory (c.f. det doc)
    for (int i = 0; i < G0_tau.size(); ++i) {
      dets.emplace_back(G0hat_t{G0_tau[i], params.alpha[i]}, params.det_init_size);
      dets.back().set_singular_threshold(params.det_singular_threshold);
      dets.back().set_n_operations_before_check(params.det_n_operations_before_check);
      dets.back().set_precision_warning(params.det_precision_warning);
      dets.back().set_precision_error(params.det_precision_error);
    }
  }

} // namespace triqs_ctint
