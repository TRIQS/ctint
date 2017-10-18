#pragma once
#include "./types.hpp"
#include "./vertex.hpp"
#include "./params.hpp"
#include "./dets.hpp"
#include <triqs/det_manip.hpp>
#include <triqs/mc_tools.hpp>

namespace triqs_ctint {

  /// Type of the Monte-Carlo configuration
  class qmc_config_t {
    public:
    qmc_config_t(params_t const &params, block_gf<imtime, matrix_valued> const &G0_tau);

    /// Unordered list of all vertices currently inserted
    std::vector<vertex_t> vertex_lst;

    /// List containing the determinant for each block
    std::vector<det_t> dets;

    /// Current perturbation order
    int perturbation_order() const { return vertex_lst.size(); }
  };

} // namespace triqs_ctint
