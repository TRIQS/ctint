#include "../types.hpp"
#include "./spinflip.hpp"

namespace triqs_ctint::moves {

  mc_weight_t spinflip::attempt() {
    TRIQS_ASSERT(n_spinflips == 1); // only single spinflip supported right now
    vpos.reserve(n_spinflips);
    vpos.clear();

    // Don't spinflip if there are not enough vertices
    if (qmc_config->vertex_lst.size() < n_spinflips) return 0.0;

    // Choose one of the vertices for spinflip
    for (int i = 0; i < n_spinflips; ++i) {
      int this_vpos = rng(qmc_config->vertex_lst.size());
      // In case of double spinflip we pick up another vertex
      if (std::find(vpos.begin(), vpos.end(), this_vpos) == vpos.end()) {
        vpos.emplace_back(this_vpos);
        // Lazy insert and capture the weight for the vertex
        auto &v = qmc_config->vertex_lst[this_vpos];
        lazy_op << v;
      } else {
        lazy_op.reset();
        return 0;
      }
    }

    // Execute the spinflip move
    g_tau_scalar_t det_ratio = lazy_op.execute_try_change_col_row();

    // Return the overall weight
    return mc_weight_t{det_ratio};
  }

  mc_weight_t spinflip::accept() {
    for (auto &d : qmc_config->dets) d.complete_operation();

    // Flip spin in the vertex list
    for (auto const &this_vpos : vpos) {
      auto &v = qmc_config->vertex_lst[this_vpos];
      v.s = 1 - v.s;
    }

    return 1.0; // no need for a correction of the sign
  }

  void spinflip::reject() {
    for (auto &d : qmc_config->dets) d.reject_last_try(); // reject the last try in all determinants
  }

} // namespace triqs_ctint::moves
