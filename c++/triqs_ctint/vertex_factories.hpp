#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <vector>
#include "./vertex.hpp"
#include "./params.hpp"
#include "./solver_core.hpp"

namespace triqs_ctint {

  /// A vertex factory that can generate random vertices for the CTInt
  using vertex_factory_t = std::function<vertex_t()>;

  /// Function that returns a list vertex factories, one for each interaction type
  std::vector<vertex_factory_t> make_vertex_factories(params_t const &params, triqs::mc_tools::random_generator &rng,
                                                      std::optional<block2_gf_const_view<imtime, matrix_valued>> D0_tau,
                                                      std::optional<gf_const_view<imtime, matrix_valued>> Jperp_tau);

} // namespace triqs_ctint
