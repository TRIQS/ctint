// Copyright (c) 2017--present, The Simons Foundation
// This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

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
                                                      std::optional<block2_gf_const_view<imfreq, matrix_valued>> D0_iw,
                                                      std::optional<gf_const_view<imfreq, matrix_valued>> Jperp_iw);

} // namespace triqs_ctint
