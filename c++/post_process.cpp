#include "post_process.hpp"
#include <triqs/gfs.hpp>
#include <triqs/gfs/singularity/fit_tail.hpp>

//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

namespace triqs_ctint {

  void Giw_from_M_iw(block_gf<imfreq, matrix_valued> &M_iw, block_gf<imfreq, matrix_valued> &G0_shift_iw, block_gf<imfreq, matrix_valued> &Giw) {
    //make Giw from g0 and M_iw
    triqs::clef::placeholder<0> w_;
    for (int b = 0; b < M_iw.size(); b++) { Giw[b](w_) << G0_shift_iw[b](w_) + G0_shift_iw[b](w_) * M_iw[b](w_) * G0_shift_iw[b](w_); }
  }

  void Sigma_iw_from_M_iw(block_gf<imfreq, matrix_valued> &M_iw, block_gf<imfreq, matrix_valued> &G0_shift_iw, block_gf<imfreq, matrix_valued> &Giw,
                          block_gf<imfreq, matrix_valued> &Sigma_iw, std::vector<std::vector<double>> const &fact) {
    for (int b = 0; b < M_iw.size(); b++) {
      triqs::clef::placeholder<0> w_;
      matrix<dcomplex> m_shift(M_iw[b].target_shape());
      m_shift() = 0;
      for (int i = 0; i < fact[b].size(); i++) m_shift(i, i) = fact[b][i];
      Sigma_iw[b](w_) << (1 / Giw[b](w_)) * M_iw[b](w_) * G0_shift_iw[b](w_) - m_shift;
    }
  }

} // namespace triqs_ctint
