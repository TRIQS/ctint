#include "./average_k.hpp"

namespace triqs_ctint::measures {

  average_k::average_k(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : average_k_(results->average_k), qmc_config(qmc_config_) {
    average_k_ = 0.0;
  }

  void average_k::accumulate(mc_weight_t sign) {
    average_k_ += qmc_config.vertex_lst.size();
    ++count;
  }

  void average_k::collect_results(triqs::mpi::communicator const &comm) {
    average_k_ = mpi_all_reduce(average_k_, comm);
    count      = mpi_all_reduce(count, comm);
    average_k_ = average_k_ / count;
  }

} // namespace triqs_ctint::measures
