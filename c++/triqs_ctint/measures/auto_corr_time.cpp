#include "./auto_corr_time.hpp"

namespace triqs_ctint::measures {

  auto_corr_time::auto_corr_time(params_t const &, qmc_config_t const &qmc_config, container_set *results)
    : qmc_config(qmc_config), auto_corr_time_(results->auto_corr_time) {}

  void auto_corr_time::accumulate(mc_weight_t sign) {
    log_accs[0] << sign;
    log_accs[1] << qmc_config.perturbation_order();
  }

  void auto_corr_time::collect_results(mpi::communicator const &comm) {
    using namespace triqs::stat;

    auto_corr_time_ = 0.0;

    for (auto &log_acc : log_accs) {
      auto [errs, counts] = log_acc.log_bin_errors_all_reduce(comm);

      // Estimate auto-correlation time
      if (comm.rank() == 0 && errs[0] > 0) {
	auto_corr_time_ = std::max(auto_corr_time_, tau_estimate_from_errors(errs[int(0.7 * errs.size())], errs[0]));
      }

      // Reset the accumulator
      log_acc = {0.0, -1, 0};
    }
    mpi::broadcast(auto_corr_time_, comm, 0);
  }

} // namespace triqs_ctint::measures
