#include <triqs_ctint/nfft_buf.hpp>
#include <random>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs::utility;

// Check that int holds at least 4 Bytes
static_assert(sizeof(int) >= sizeof(int32_t), " Error: sizeof(int) < 4 Byte ");

/********************* Fixture Common to all TEST_F ********************/
class Nfft : public ::testing::Test {

  protected:
  // parameters
  double beta = 20.0;
  int n_iw    = 100;

  // gf containers common to multiple tests
  std::array<long, 2> shape;
  gf<imfreq, matrix_valued> giw_exact;

  void SetUp() override {
    shape = std::array<long, 2>{1, 1};

    // Init exact reference gf
    nda::clef::placeholder<0> iw_;
    giw_exact = gf<imfreq, matrix_valued>{{beta, Fermion, n_iw}, shape};
    giw_exact(iw_) << 1.0 / (iw_ - 1.0);
  }

  // function to be transformed
  double f_tau(double tau) { return -std::exp(-tau) / (1 + std::exp(-beta)); } // bosonic: std::exp(-tau) (std::exp(-beta) - 1)
};

/********************* EQUIDISTANT TRANSFORM ********************/
TEST_F(Nfft, Equid) { // NOLINT

  // Parameters
  int n_tau    = 10000;
  int buf_size = n_tau;

  // Create container for Gf from equid nfft
  auto giw_nfft_equid = gf<imfreq, matrix_valued>{{beta, Fermion, n_iw}, shape};

  // nfft_buffer
  nfft_buf_t<1> buf_equid(giw_nfft_equid.data()(range::all, 0, 0), buf_size, beta);

  // Generate data with equidistant \tau-grid (care for weights at beginning and end)
  buf_equid.push_back({0.0}, 0.5 * f_tau(0.0));
  for (int i = 1; i < n_tau - 1; ++i) {
    double tau = beta * i / (n_tau - 1);
    buf_equid.push_back({tau}, f_tau(tau));
  }
  buf_equid.push_back({beta - 1e-10}, 0.5 * f_tau(beta - 1e-10));

  // normalize, and care for half-points at 0^+ and \beta^-
  giw_nfft_equid *= beta / (n_tau - 1);

  // Generate Gf with fftw
  auto gtau    = gf<imtime, matrix_valued>{{beta, Fermion, n_tau}, shape};
  auto iw_mesh = make_adjoint_mesh(gtau.mesh(), n_iw);
  for (auto tau : gtau.mesh()) gtau[tau] = f_tau(tau);
  auto giw_fftw = make_gf_from_fourier(gtau, iw_mesh, make_zero_tail(gtau));

  // Compare to exact and fftw
  EXPECT_GF_NEAR(giw_nfft_equid, giw_exact, 1e-4); // Same order of fftw discretization error
  EXPECT_GF_NEAR(giw_nfft_equid, giw_fftw, 1e-8);  // Only small deviation due to truncation/oversampling factor (see Fig.3 Notes Josef)

  // -- Now check multiple transforms

  // Create container for Gf from multi nfft
  auto giw_nfft_multi = gf<imfreq, matrix_valued>{{beta, Fermion, n_iw}, shape};

  // nfft_buffer with size 4000 = 2 * buf_size / 5
  nfft_buf_t<1> buf_multi(giw_nfft_multi.data()(range::all, 0, 0), 2 * buf_size / 5, beta);

  // Generate data with equidistant \tau-grid (care for weights at beginning and end)
  // Buffer performs multiple transforms as buf_size < n_tau
  buf_multi.push_back({0.0}, 0.5 * f_tau(0.0));
  for (int i = 1; i < n_tau - 1; ++i) {
    double tau = beta * i / (n_tau - 1);
    buf_multi.push_back({tau}, f_tau(tau));
  }
  buf_multi.push_back({beta - 1e-10}, 0.5 * f_tau(beta - 1e-10));

  // Care to flush remaining points since 10000 is not an even multiple of 4000
  buf_multi.flush();

  // normalize, and care for half-points at 0^+ and \beta^-
  giw_nfft_multi *= beta / (n_tau - 1);

  // Compare one-shot with multi-nfft
  EXPECT_GF_NEAR(giw_nfft_multi, giw_nfft_equid, 1e-12); // There should be no difference between the one-shot and multi nfft

  // Write to file
  // h5::file h5file("nfft.equid.h5", 'w');
  // h5_write(h5file "giw_nfft_equid", giw_nfft_equid);
  // h5_write(h5file "giw_fftw", giw_fftw);
  // h5_write(h5file "giw_exact", giw_exact);
  // h5_write(h5file "giw_nfft_multi", giw_nfft_multi);
}

/*********************  NON-EQUIDISTANT TRANSFORM ********************/
TEST_F(Nfft, Rng) { // NOLINT

  // Parameters
  int n_tau    = 1e+6;
  int buf_size = n_tau;

  // std random generator
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  // Create container for Gf from nfft
  auto giw_nfft_rng = gf<imfreq, matrix_valued>{{beta, Fermion, n_iw}, shape};

  // nfft_buffer
  n_tau    = 1e+6;
  buf_size = n_tau;
  nfft_buf_t<1> buf_rng(giw_nfft_rng.data()(range::all, 0, 0), buf_size, beta);

  // Generate values at random tau points
  for (int i = 0; i < n_tau; ++i) {
    double tau = double(distribution(generator)) * beta;
    buf_rng.push_back({tau}, f_tau(tau));
  }

  // normalize
  giw_nfft_rng *= beta / n_tau;

  // Compare
  EXPECT_GF_NEAR(giw_nfft_rng, giw_exact, 1e-2); // Expect a Monte-Carlo Error of order 1/sqrt(n_tau)

  // Write to file
  // h5::file h5file("nfft.rng.h5", 'w');
  // h5_write(h5file "giw_nfft_rng", giw_nfft_rng);
  // h5_write(h5file "giw_exact", giw_exact);
}

/*********************  2D NFFT TRANSFORM ********************/
TEST_F(Nfft, 2D) { // NOLINT

  // Parameters
  int n_tau    = 6 * n_iw + 1;
  int buf_size = n_tau * n_tau;

  // Create container for Gf from nfft
  auto giw_nfft_2d = gf<prod<imfreq, imfreq>>{{{beta, Fermion, n_iw}, {beta, Fermion, n_iw}}, shape};

  // nfft_buffer
  nfft_buf_t<2> buf_2d(slice_target_to_scalar(giw_nfft_2d, 0, 0).data(), buf_size, beta);

  // ==== Generate 2d data with equidistant \tau-grids (care for weights at edges and corners)
  // Corner Points with weight 0.25
  buf_2d.push_back({0.0, 0.0}, 0.25 * f_tau(0.0) * f_tau(0.0));
  buf_2d.push_back({0.0, beta - 1e-10}, 0.25 * f_tau(0.0) * f_tau(beta - 1e-10));
  buf_2d.push_back({beta - 1e-10, 0.0}, 0.25 * f_tau(beta - 1e-10) * f_tau(0.0));
  buf_2d.push_back({beta - 1e-10, beta - 1e-10}, 0.25 * f_tau(beta - 1e-10) * f_tau(beta - 1e-10));
  for (int i = 1; i < n_tau - 1; ++i) {
    double tau_i = beta * i / (n_tau - 1);
    for (int j = 1; j < n_tau - 1; ++j) {
      double tau_j = beta * j / (n_tau - 1);
      // Core points with full weight
      buf_2d.push_back({tau_i, tau_j}, f_tau(tau_i) * f_tau(tau_j));
    }
    // Edge points with weight 0.5
    buf_2d.push_back({tau_i, 0.0}, 0.5 * f_tau(tau_i) * f_tau(0.0));
    buf_2d.push_back({tau_i, beta - 1e-10}, 0.5 * f_tau(tau_i) * f_tau(beta - 1e-10));
  }
  for (int j = 1; j < n_tau - 1; ++j) {
    double tau_j = beta * j / (n_tau - 1);
    // Edge points with weight 0.5
    buf_2d.push_back({0, tau_j}, 0.5 * f_tau(0) * f_tau(tau_j));
    buf_2d.push_back({beta - 1e-10, tau_j}, 0.5 * f_tau(beta - 1e-10) * f_tau(tau_j));
  }
  // ====

  // normalize, care for half-points at 0^+ and \beta^-
  giw_nfft_2d *= beta * beta / (n_tau - 1) / (n_tau - 1);

  // === Generate 2d Gf with fftw
  // Create 1d giw from fftw
  auto gtau2   = gf<imtime, matrix_valued>{{beta, Fermion, n_tau}, shape};
  auto iw_mesh = make_adjoint_mesh(gtau2.mesh(), n_iw);
  for (auto tau : gtau2.mesh()) gtau2[tau] = f_tau(tau);
  auto giw_fftw = make_gf_from_fourier(gtau2, iw_mesh, make_zero_tail(gtau2));
  // Create giw_fftw_2d from product of giw_fftw
  auto giw_fftw_2d = gf<prod<imfreq, imfreq>>{{{beta, Fermion, n_iw}, {beta, Fermion, n_iw}}, shape};
  for (auto iw1 : giw_fftw.mesh())
    for (auto iw2 : giw_fftw.mesh()) giw_fftw_2d[iw1, iw2] = giw_fftw[iw1] * giw_fftw[iw2];

  // Init exact reference gf
  nda::clef::placeholder<0> iw1_;
  nda::clef::placeholder<1> iw2_;
  auto giw_exact_2d = gf<prod<imfreq, imfreq>>{{{beta, Fermion, n_iw}, {beta, Fermion, n_iw}}, shape};
  giw_exact_2d(iw1_, iw2_) << 1.0 / (iw1_ - 1.0) / (iw2_ - 1.0);

  // Compare
  EXPECT_GF_NEAR(giw_nfft_2d, giw_exact_2d, 1e-2); // Same order of fftw discretization error
  EXPECT_GF_NEAR(giw_nfft_2d, giw_fftw_2d, 1e-8);  // Only small deviation due to truncation/oversampling factor (see Fig.3 Notes Josef)

  // Write to file
  // h5::file h5file("nfft.2d.h5", 'w');
  // h5_write(h5file "arr_nfft_2d", giw_nfft_2d.data());
  // h5_write(h5file "arr_fftw_2d", giw_fftw_2d.data());
  // h5_write(h5file "arr_exact_2d", giw_exact_2d.data());
}

MAKE_MAIN;
