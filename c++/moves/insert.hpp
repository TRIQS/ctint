#pragma once
#include "./qmc_config.hpp"
#include "./lazy_det_operation.hpp"
#include "./vertex_factories.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace triqs_ctint::moves {

  /// The move that inserts one or multiple vertices into the configuration
  struct insert {

    /// The Monte-Carlo configuration
    qmc_config_t *qmc_config;

    /// Factory that randomly generates vertices
    std::vector<vertex_factory_t> const &vertex_factories;

    /// The random number generator
    triqs::mc_tools::random_generator &rng;

    /// Switch for double vertex insertions
    bool double_insert;

    /// Object that allows to delay determinant operations, necessary for multi-inserts/removes
    lazy_det_operation_t lazy_op = lazy_det_operation_t{&qmc_config->dets};

    /// Attempt vertex insertion
    mc_weight_t attempt();

    /// Accept vertex insertion
    mc_weight_t accept();

    /// Reject vertex insertion
    void reject();
  };

} // namespace triqs_ctint::moves
