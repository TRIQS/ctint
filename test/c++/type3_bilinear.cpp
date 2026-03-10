#include <triqs_ctint/nfft_buf.hpp>
#include <triqs/mesh/matsubara_freq.hpp>
#include <random>
#include <gtest/gtest.h>

using namespace triqs::utility;

using dcomplex = std::complex<double>;

// Simple Gauss-Jordan inversion for small complex matrices used only in tests
static std::vector<std::vector<dcomplex>> invert_matrix(std::vector<std::vector<dcomplex>> A) {
  const size_t N = A.size();
  std::vector<std::vector<dcomplex>> I(N, std::vector<dcomplex>(N, dcomplex{0, 0}));
  for (size_t i = 0; i < N; ++i) I[i][i] = dcomplex{1, 0};

  for (size_t i = 0; i < N; ++i) {
    // find pivot
    size_t piv = i;
    double maxabs = std::abs(A[i][i]);
    for (size_t r = i + 1; r < N; ++r) {
      double a = std::abs(A[r][i]);
      if (a > maxabs) { maxabs = a; piv = r; }
    }
    if (piv != i) { std::swap(A[i], A[piv]); std::swap(I[i], I[piv]); }

    // normalize row
    dcomplex diag = A[i][i];
    if (std::abs(diag) == 0.0) throw std::runtime_error("Singular matrix in test invert");
    for (size_t c = 0; c < N; ++c) { A[i][c] /= diag; I[i][c] /= diag; }

    // eliminate other rows
    for (size_t r = 0; r < N; ++r) {
      if (r == i) continue;
      dcomplex fac = A[r][i];
      for (size_t c = 0; c < N; ++c) { A[r][c] -= fac * A[i][c]; I[r][c] -= fac * I[i][c]; }
    }
  }
  return I;
}

TEST(Type3Bilinear, FinufftEqualsBilinear) {
  const int N = 12; // small test size
  const double beta = 7.5;

  std::default_random_engine gen(12345);
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  std::uniform_real_distribution<double> rd(-1.0, 1.0);

  // Build a random (but well-conditioned) complex matrix G
  std::vector<std::vector<dcomplex>> G(N, std::vector<dcomplex>(N));
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) G[i][j] = dcomplex(rd(gen), rd(gen));
  // make diagonally dominant for stability
  for (int i = 0; i < N; ++i) G[i][i] += dcomplex(N, 0);

  auto Ginv = invert_matrix(G);

  // random taus and channel labels
  std::vector<double> taus(N);
  std::vector<int> c(N), cbar(N);
  for (int i = 0; i < N; ++i) {
    taus[i] = ud(gen) * beta;
    c[i] = (i % 3);      // a few channels
    cbar[i] = ((i + 1) % 3);
  }

  const int a_channel = 1;
  const int b_channel = 2;
  const int n_a = 3;
  const int n_b = -2;
  const double omega_a = (2 * n_a + 1) * M_PI / beta;
  const double omega_b = (2 * n_b + 1) * M_PI / beta;

  // build v_i = delta_{c_i,a} * e^{-i omega_a tau_i}
  std::vector<dcomplex> v(N, dcomplex{0, 0});
  for (int i = 0; i < N; ++i)
    if (c[i] == a_channel) v[i] = std::exp(dcomplex(0, -omega_a * taus[i]));

  // compute x = G^{-1} * v
  std::vector<dcomplex> x(N, dcomplex{0, 0});
  for (int j = 0; j < N; ++j)
    for (int i = 0; i < N; ++i) x[j] += Ginv[j][i] * v[i];

  // bilinear result M = sum_j delta_{cbar_j,b} * x_j * exp(i omega_b tau_j)
  dcomplex M_ref{0, 0};
  for (int j = 0; j < N; ++j)
    if (cbar[j] == b_channel) M_ref += x[j] * std::exp(dcomplex(0, omega_b * taus[j]));

  // Now compute same quantity via FINUFFT type-3 by pushing ftau = x_j * delta_{cbar_j,b}
  std::vector<triqs::mesh::matsubara_freq> target_mf;
  target_mf.emplace_back(triqs::mesh::matsubara_freq(n_b, beta, triqs::mesh::Fermion));

  nda::vector<dcomplex> fiw_out(1);
  fiw_out = 0;
  nfft_buf_t<1> buf(fiw_out, target_mf, N, nfft_type_t::type3);

  for (int j = 0; j < N; ++j) {
    dcomplex ftau = (cbar[j] == b_channel) ? x[j] : dcomplex{0, 0};
    buf.push_back({taus[j]}, ftau);
  }
  buf.flush();

  dcomplex M_finufft = fiw_out(0);

  double err = std::abs(M_finufft - M_ref);
  EXPECT_LT(err, 1e-10) << "Mismatch: finufft=" << M_finufft << " ref=" << M_ref << " err=" << err;
}
