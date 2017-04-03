#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/gfs.hpp>
#include <vector>
#include "./vertex.hpp"
#include "./params.hpp"
#include "./solver_core.hpp"

namespace triqs_ctint {

  /// A vertex factory that can generate random vertices for the CTInt
  using vertex_factory_t = std::function<vertex_t()>;

  /// Function that returns a list vertex factories, one for each interaction type
  std::vector<vertex_factory_t> make_vertex_factories(params_t const &params, triqs::mc_tools::random_generator &rng,
                                                      std::optional<block_gf_const_view<imfreq, matrix_valued>> D0_iw,
                                                      std::optional<gf_const_view<imfreq, matrix_valued>> Jperp_iw);

} // namespace triqs_ctint
