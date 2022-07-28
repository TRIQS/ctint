#include "./average_sign.hpp"

namespace triqs_ctint::measures {

  average_sign::average_sign(params_t const &, qmc_config_t const &qmc_config_, container_set *results)
     : average_sign_(results->average_sign), nmeasures(results->nmeasures), qmc_config(qmc_config_) {
    average_sign_ = 0.0;
    nmeasures     = 0;
  }

  void average_sign::accumulate(mc_weight_t sign) {
    average_sign_ += sign;

    // mc_sign_t det_static  = 0.0;
    //mc_sign_t det_dynamic = 0.0;

      //if (params.sign_analysis) {
      // - Get matrix for each block (qmc_config.dets[bl].matrix();)
      // - Set elements for equal / non-equal times to zero in the matrix
      // - Calculate det
      // - Multiply dets and take the sign
      //}
    ++count;
  }

  /**
* add sign_analysis to params -> add only sign measure
* initialize additional dets in config for static/dynamic sign
* Update insert / remove to update also other dets
*/

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
