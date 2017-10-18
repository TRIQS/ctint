#include "./average_sign.hpp"

namespace triqs_ctint::measures {

  average_sign::average_sign(params_t const &params_, qmc_config_t const &qmc_config_, container_set *results)
     : average_sign_(results->average_sign) {
    average_sign_ = 0.0;
  }

  void average_sign::accumulate(double sign) {
    average_sign_ += sign;
    ++count;
  }

  void average_sign::collect_results(triqs::mpi::communicator const &comm) {
    average_sign_ = mpi_all_reduce(average_sign_, comm);
    count         = mpi_all_reduce(count, comm);
    average_sign_ = average_sign_ / count;
  }

} // namespace triqs_ctint::measures
