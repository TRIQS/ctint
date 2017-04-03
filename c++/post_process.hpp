#pragma once
#include "./qmc_config.hpp"
#include "./params.hpp"

namespace triqs_ctint {

  void Giw_from_M_iw(block_gf<imfreq, matrix_valued> &M_iw, block_gf<imfreq, matrix_valued> &G0_shift_iw, block_gf<imfreq, matrix_valued> &Giw);
  void Sigma_iw_from_M_iw(block_gf<imfreq, matrix_valued> &M_iw, block_gf<imfreq, matrix_valued> &G0_shift_iw, block_gf<imfreq, matrix_valued> &Giw,
                          block_gf<imfreq, matrix_valued> &Sigmaw, std::vector<std::vector<double>> const &fact);

} // namespace triqs_ctint
