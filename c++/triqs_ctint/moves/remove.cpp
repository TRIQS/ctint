#include "./remove.hpp"

namespace triqs_ctint::moves {

  mc_weight_t remove::attempt() {

    // Don't remove if config is empty
    if (qmc_config->perturbation_order() == 0) return 0.0;

    // Prepare remove of single vertex
    auto single_remove = [&](int vpos) {

      // Lazy remove the vertex at given position
      auto &v = qmc_config->vertex_lst[vpos];
      lazy_op << v;

      // Calculate insertion proposition probability
      double insert_proposition_proba = v.proposition_proba / vertex_factories.size();

      // Return product of vertex amplitude and insertion proposition probability
      return 1.0 / v.amplitude * insert_proposition_proba;
    };

    // Choose one of the vertices for removal
    vpos = rng(qmc_config->perturbation_order());

    // Lazy insert and capture the weight for the first vertex
    U_scalar_t ratio = single_remove(vpos);

    // In case of double remove we pick up another vertex
    if (double_remove) {
      vpos2 = rng(qmc_config->perturbation_order());
      if (vpos == vpos2) return 0;
      ratio *= single_remove(vpos2);
    }

    // Execute the removal move
    g_tau_scalar_t det_ratio = lazy_op.execute_try_remove();

    // Calculate the removal proposition probability
    double remove_proposition_proba = 1.0 / (qmc_config->perturbation_order() * (double_remove ? qmc_config->perturbation_order() : 1));

    // Return the overall weight
    return mc_weight_t{det_ratio} * ratio / remove_proposition_proba;
  }

  mc_weight_t remove::accept() {

    // Finish the removal from the determinants
    for (auto &d : qmc_config->dets) d.complete_operation(); // maybe doing nothing if not try_insert but it is quick

    // Remove the vertices from vertex_lst
    if (double_remove) {
      if (vpos > vpos2) std::swap(vpos, vpos2); // ensure vpos < vpos2
      qmc_config->vertex_lst.erase(begin(qmc_config->vertex_lst) + vpos2);
    }
    qmc_config->vertex_lst.erase(begin(qmc_config->vertex_lst) + vpos);

    return 1.0; // no need for a correction of the sign
  }

  void remove::reject() {
    for (auto &d : qmc_config->dets) d.reject_last_try(); // reject the last try in all determinants
  }

} // namespace triqs_ctint::moves
