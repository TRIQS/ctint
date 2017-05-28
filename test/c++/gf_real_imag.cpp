#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>
#include "types.hpp"

const std::complex<double> I(0.0, 1.0);

using namespace triqs::gfs;
using namespace triqs::clef;

/********************* EQUIDISTANT TRANSFORM ********************/
TEST(Gfs, Real_Imag) {

  // Parameters
  int n_iw    = 100;
  double beta = 10.0; 

  placeholder<0> iw_; 

  auto G = gf<imfreq, matrix_valued>{{beta, Fermion, n_iw}, make_shape(2,2)}; 

  G[iw_] << 1.0 / ( iw_ + 4.0 ); 

  //gf<imfreq, matrix_valued> G_Res = get_real(G) + I * get_imag(G);  // FIXME OPERATOR= not implemented

  //EXPECT_GF_NEAR(G, G_Res, 1e-15);  // Only small deviation due to truncation/oversampling factor (see Fig.3 Notes Josef)
}

MAKE_MAIN;
