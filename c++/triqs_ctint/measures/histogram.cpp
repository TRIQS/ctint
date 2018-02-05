#include "./histogram.hpp"

namespace triqs_ctint::measures {

  histogram::histogram(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : histogram_(results->histogram), qmc_config(qmc_config_) {
    results->histogram.resize(1000);
  }

  void histogram::accumulate(mc_weight_t sign) {
    int k = qmc_config.perturbation_order();
    if (k >= histogram_.size()) TRIQS_RUNTIME_ERROR << "Expansion order larger than " << histogram_.size() << ": exiting.";
    histogram_[k] += 1.;
    ++N;
  }

  void histogram::collect_results(triqs::mpi::communicator const &comm) {
    N = mpi_all_reduce(N, comm);
    for (auto &h_k : histogram_) {
      h_k = mpi_all_reduce(h_k, comm);
      h_k = h_k / N;
    }
  }
} // namespace triqs_ctint::measures
