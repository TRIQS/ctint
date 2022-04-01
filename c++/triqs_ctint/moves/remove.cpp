#include "./remove.hpp"

namespace triqs_ctint::moves {

  mc_weight_t remove::attempt() {
    TRIQS_ASSERT(n_removals > 0);
    vpos.reserve(n_removals);
    vpos.clear();

    // Don't remove if config is empty
    if (qmc_config->perturbation_order() == 0) return 0.0;

    // Prepare remove of single vertex
    auto single_remove = [&](int p) {
      // Lazy remove the vertex at given position
      auto &v = qmc_config->vertex_lst[p];
      lazy_op << v;

      // Calculate insertion proposition probability
      double insert_proposition_proba = v.proposition_proba / vertex_factories.size();

      // Return ratio of insertion proposition probability and vertex amplitude
      return insert_proposition_proba / v.amplitude;
    };

    // Choose one of the vertices for removal
    U_scalar_t ratio = 1;
    for (int i = 0; i < n_removals; ++i) {
      int const this_vpos = rng(qmc_config->perturbation_order());
      // In case of double remove we pick up another vertex
      if (std::find(vpos.begin(), vpos.end(), this_vpos) == vpos.end()) {
        vpos.emplace_back(this_vpos);
        // Lazy insert and capture the weight for the vertex
        ratio *= single_remove(this_vpos);
      } else {
        lazy_op.reset();
        return 0;
      }
    }

    // Execute the removal move
    g_tau_scalar_t det_ratio = lazy_op.execute_try_remove();


    double remove_proposition_proba = 0.0;
    if(max_order == -1 || qmc_config->perturbation_order() < max_order) {
      remove_proposition_proba = 1.0 / std::pow(qmc_config->perturbation_order(), n_removals);
    }

    // Return the overall weight
    return mc_weight_t{det_ratio} * ratio / remove_proposition_proba;
  }

  mc_weight_t remove::accept() {

    // Finish the removal from the determinants
    for (auto &d : qmc_config->dets) d.complete_operation(); // maybe doing nothing if not try_insert but it is quick

    // sort the vertices in reverse order
    std::sort(vpos.begin(), vpos.end(), std::greater<>());

    // Remove the vertices from vertex_lst
    for (auto &&this_vpos : vpos) {
      qmc_config->vertex_lst.erase(begin(qmc_config->vertex_lst) + this_vpos);
    }

    return 1.0; // no need for a correction of the sign
  }

  void remove::reject() {
    for (auto &d : qmc_config->dets) d.reject_last_try(); // reject the last try in all determinants
  }

} // namespace triqs_ctint::moves
