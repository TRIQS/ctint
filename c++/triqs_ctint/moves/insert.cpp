#include "./insert.hpp"

namespace triqs_ctint::moves {

  mc_weight_t insert::attempt() {

    auto single_insert = [&]() {

      // Pick up one factory at random
      int n_fact       = vertex_factories.size();
      auto const &fact = vertex_factories[rng(n_fact)];

      // Generate a vertex from the factory
      vertex_t v = fact();

      // Lazy add the vertex to the det
      lazy_op << v;

      // Calculate Monte-Carlo weight
      double insert_proposition_proba = v.proposition_proba / n_fact;
      U_scalar_t weight               = v.amplitude / insert_proposition_proba;

      // Move new vertex to its storage
      qmc_config->vertex_lst.push_back(std::move(v));

      return weight;
    };

    // Lazy insert and capture the weight for one vertex
    U_scalar_t ratio = single_insert();

    // Perform previous step again in case of double removal
    if (double_insert) ratio *= single_insert();

    // Execute the insertion move
    g_tau_scalar_t det_ratio = lazy_op.execute_try_insert();

    // Calculate the removal proposition probability
    double remove_proposition_proba = 1.0 / (qmc_config->perturbation_order() * (double_insert ? qmc_config->perturbation_order() : 1));

    // Return the overall weight
    return mc_weight_t{det_ratio} * remove_proposition_proba * ratio;
  }

  mc_weight_t insert::accept() {
    for (auto &d : qmc_config->dets) d.complete_operation(); // maybe doing nothing if not try_insert but it is quick
    return 1.0;                                              // no need for a correction of the sign
  }

  void insert::reject() {
    for (int i = 0; i < (double_insert ? 2 : 1); ++i) qmc_config->vertex_lst.pop_back(); // remove the last insertions.
    for (auto &d : qmc_config->dets) d.reject_last_try();                                // reject the last try in all determinants
  }

} // namespace triqs_ctint::moves
