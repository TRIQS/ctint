#pragma once
#include "./qmc_config.hpp"
#include "./lazy_det_operation.hpp"
#include "./vertex_factories.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace triqs_ctint::moves {

  /// The move that removes one or multiple vertices from the configuration
  struct remove {

    /// The Monte-Carlo configuration
    qmc_config_t *qmc_config;

    /// Factory that randomly generates vertices
    std::vector<vertex_factory_t> const &vertex_factories;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// Switch for double vertex removals
    double double_remove;

    int vpos = 0, vpos2 = 0;

    /// Object that allows to delay determinant operations, necessary for multi-inserts/removes
    lazy_det_operation_t lazy_op = lazy_det_operation_t{&qmc_config->dets};

    /// Attempt vertex removal
    mc_weight_t attempt();

    /// Accept vertex removal
    mc_weight_t accept();

    /// Reject vertex removal
    void reject();
  };

} // namespace triqs_ctint::moves
