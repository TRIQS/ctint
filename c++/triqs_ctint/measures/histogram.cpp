#include "./histogram.hpp"

namespace triqs_ctint::measures {

  histogram::histogram(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : qmc_config(qmc_config_), histogram_(results->histogram) {
    results->histogram = std::vector<double>(4);
  }

  void histogram::accumulate(mc_weight_t sign) {
    int k = qmc_config.perturbation_order();
    while (k >= histogram_->size()) histogram_->resize(2 * histogram_->size());
    histogram_.value()[k] += 1.;
    ++N;
  }

  void histogram::collect_results(mpi::communicator const &comm) {
    N = mpi::all_reduce(N, comm);

    // Make sure that all mpi threads have an equally sized histogram
    auto max_k_vec         = std::vector<size_t>(comm.size());
    max_k_vec[comm.rank()] = histogram_->size();
    max_k_vec              = mpi::all_reduce(max_k_vec, comm);
    histogram_->resize(*std::max_element(max_k_vec.begin(), max_k_vec.end()));

    // Reduce histogram over all mpi threads
    histogram_ = mpi::all_reduce(histogram_.value(), comm);
    for (auto &h_k : histogram_.value()) h_k = h_k / N;
  }
} // namespace triqs_ctint::measures
