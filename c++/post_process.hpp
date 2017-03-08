#pragma once
#include "./qmc_config.hpp"
#include "./params.hpp"

namespace triqs_ctint {

  /// Calculate the vertex function $F$ from the the building blocks $M4_iw$ and $M_iw$
  chi4_iw_t F_from_M4(chi4_iw_t::const_view_type M4_iw, block_gf_const_view<imfreq, matrix_valued> M_iw, block_gf_const_view<imfreq, matrix_valued> G0_iw);

  /// Calculate the $K_2$ function from the building block $M3_iw$ and $M_iw$
  chi3_iw_t K2_from_M3(chi3_iw_t::const_view_type M3_iw, g_iw_t::const_view_type M_iw, g_iw_t::const_view_type G0_iw); 

} // namespace triqs_ctint
