#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs_ctint/types.hpp>


using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace triqs::clef;

TEST(Gfs, Real_Imag) { // NOLINT

  // Parameters
  int n_iw    = 100;
  double beta = 10.0; 

  placeholder<0> iw_; 

  auto G = gf<imfreq, matrix_valued>{{beta, Fermion, n_iw}, make_shape(2,2)}; 

  G[iw_] << 1.0 / ( iw_ + 4.0 ); 
  
  //auto G_real = real(G); 
  //auto G_imag = imag(G); 

  // FIXME 
  //const std::complex<double> I(0.0, 1.0);
  //auto G_Res(iw_) <<  I * G_imag(iw_) + G_real(iw_);

  //EXPECT_GF_NEAR(G, G_Res);
}

MAKE_MAIN;
