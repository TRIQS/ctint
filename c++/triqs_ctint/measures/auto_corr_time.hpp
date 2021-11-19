#pragma once
#include <triqs/stat/accumulator.hpp>
#include "../qmc_config.hpp"
#include "../container_set.hpp"

namespace triqs_ctint::measures {

  /// Measurement auto-correlation time based on the partition function
  struct auto_corr_time {

    auto_corr_time(params_t const &, qmc_config_t const &, container_set *results);

    /// Accumulate average sign
    void accumulate(mc_weight_t sign);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &comm);

    private:
    // Reference to double for accumulation
    double &auto_corr_time_;

    // Accumulator
    triqs::stat::accumulator<mc_weight_t> log_acc = {0.0, -1, 0};
  };

} // namespace triqs_ctint::measures
