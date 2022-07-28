#include "./average_k.hpp"

namespace triqs_ctint::measures {

  average_k::average_k(params_t const &, qmc_config_t const &qmc_config_, container_set *results)
     : qmc_config(qmc_config_), average_k_(results->average_k) {
    average_k_ = 0.0;
  }

  void average_k::accumulate(mc_weight_t) {
    average_k_ += qmc_config.perturbation_order();
    ++N;
  }

  void average_k::collect_results(mpi::communicator const &comm) {
    average_k_ = mpi::all_reduce(average_k_, comm);
    N          = mpi::all_reduce(N, comm);
    average_k_ = average_k_ / N;
  }

} // namespace triqs_ctint::measures
