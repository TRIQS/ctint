#include "./average_sign.hpp"

namespace triqs_ctint::measures {

  average_sign::average_sign(params_t const &, qmc_config_t const &, container_set *results)
     : average_sign_(results->average_sign), nmeasures(results->nmeasures) {
    average_sign_ = 0.0;
    nmeasures     = 0;
  }

  void average_sign::accumulate(mc_weight_t sign) {
    average_sign_ += sign;
    ++count;
  }

  void average_sign::collect_results(mpi::communicator const &comm) {
    average_sign_ = mpi::all_reduce(average_sign_, comm);
    count         = mpi::all_reduce(count, comm);
    average_sign_ = average_sign_ / count;
    nmeasures     = count;
  }

  std::string average_sign::report() const {
    std::ostringstream os;
    os << "Average sign: " << average_sign_ / count;
    return os.str();
  }

} // namespace triqs_ctint::measures
