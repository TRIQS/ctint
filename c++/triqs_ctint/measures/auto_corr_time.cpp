#include "./auto_corr_time.hpp"

namespace triqs_ctint::measures {

  auto_corr_time::auto_corr_time(params_t const &, qmc_config_t const &, container_set *results) : auto_corr_time_(results->auto_corr_time) {}

  void auto_corr_time::accumulate(mc_weight_t sign) { log_acc << sign; }

  void auto_corr_time::collect_results(mpi::communicator const &comm) {
    using namespace triqs::stat;

    auto [errs, counts] = log_acc.log_bin_errors_all_reduce(comm);

    // Estimate auto-correlation time
    auto_corr_time_ = 0.0;
    if (comm.rank() == 0 && errs[0] > 0) auto_corr_time_ = std::max(0.0, tau_estimate_from_errors(errs[int(0.7 * errs.size())], errs[0]));
    mpi::broadcast(auto_corr_time_, comm, 0);

    // Reset the accumulator
    log_acc = {0.0, -1, 0};
  }

} // namespace triqs_ctint::measures
