#include <triqs_ctint/solver_core.hpp>

#include <h5/h5.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs_ctint;

TEST(CtInt, Plaquette) { // NOLINT
  // --------- physical parameters ----------
  double const U    = 1.0;     // Density-density interaction
  double const t    = 1.0;     // Hopping
  double const mu   = U / 2.0; // Chemical Potential
  double const beta = 100;     // Inverse temperature

  // --------- simulation parameters ----------
  int const n_cyc = 50;

  // --------- Define hopping matrix and interaction hamiltonian ----------

  long const Nx    = 2;
  long const Ny    = 2;
  long const n_orb = Nx * Ny;

  auto hloc0_mat = nda::matrix<double>::zeros(n_orb, n_orb);

  // diagonal terms
  diagonal(hloc0_mat) -= mu;

  auto hloc0_arr = reshape(hloc0_mat, Nx, Ny, Nx, Ny);
  auto _         = nda::range::all;
  for (long x : range(0, Nx - 1)) {
    diagonal(hloc0_arr(x, _, x + 1, _)) -= t;
    diagonal(hloc0_arr(x + 1, _, x, _)) -= t;
  }
  for (long y : range(0, Ny - 1)) {
    diagonal(hloc0_arr(_, y, _, y + 1)) -= t;
    diagonal(hloc0_arr(_, y + 1, _, y)) -= t;
  }

  many_body_operator h_int;
  for (long i = 0; i < n_orb; ++i) { h_int += U * n("up", i) * n("dn", i); }

  // --------- set up block structure ---------

  gf_struct_t gf_struct{{"dn", n_orb}, {"up", n_orb}};

  // --------- Construct the ctint solver ----------
  constr_params_t pc;
  pc.beta      = beta;
  pc.gf_struct = gf_struct;
  pc.n_iw      = 100;
  pc.n_tau     = 201;

  solver_core S(pc);

  // --------- Initialize the non-interacting Green's function ----------
  for (auto &&g_bl : S.G0_iw) {
    g_bl(iw_) << iw_ - hloc0_mat;
    g_bl = inverse(g_bl);
  }

  // Define the alpha tensor
  long n_terms = std::distance(h_int.begin(), h_int.end());
  alpha_t alpha(n_terms, 2, 2, 2);
  double const delta = 0.1;
  for (long l = 0; l < n_terms; ++l) {
    alpha(l, _, _, 0) = nda::matrix<double>{{0.5 - delta, 0.}, {0., 0.5 + delta}};
    alpha(l, _, _, 1) = nda::matrix<double>{{0.5 + delta, 0.}, {0., 0.5 - delta}};
  };

  // --------- Solve! ----------

  solve_params_t ps;
  ps.h_int              = h_int;
  ps.n_s                = 2;
  ps.alpha              = alpha;
  ps.n_cycles           = n_cyc;
  ps.length_cycle       = 50;
  ps.n_warmup_cycles    = 100;
  ps.random_seed        = 34788;
  ps.measure_histogram  = true;
  ps.measure_density    = true;
  ps.measure_M_tau      = true;
  ps.measure_M_iw       = true;
  ps.measure_M4_iw      = true;
  ps.measure_M4pp_iw    = true;
  ps.measure_M4ph_iw    = true;
  ps.n_iw_M4            = 2;
  ps.n_iW_M4            = 2;
  ps.nfft_buf_size      = 100000;
  ps.measure_M3pp_iw    = true;
  ps.measure_M3ph_iw    = true;
  ps.measure_M3pp_tau   = true;
  ps.measure_M3ph_tau   = true;
  ps.measure_M3xph_tau  = true;
  ps.n_iw_M3            = 4;
  ps.n_iW_M3            = 4;
  ps.n_tau_M3           = 4;
  ps.measure_chi2pp_tau = true;
  ps.measure_chi2ph_tau = true;
  ps.n_iw_chi2          = 10;
  ps.n_tau_chi2         = 21;
  ps.measure_chiAB_tau  = true;
  ps.chi_A_vec          = {n("up", 0) + n("dn", 0)};
  ps.chi_B_vec          = {n("up", 0) + n("dn", 0)};
  ps.post_process       = false;

  S.solve(ps);

  // -------- Save in archive ---------
  auto archive = h5::file("plaquette.out.h5", 'w');
  h5_write(archive, "histogram", S.histogram);
  h5_write(archive, "density", S.density);
  h5_write(archive, "M_tau", S.M_tau);
  h5_write(archive, "M_iw", S.M_iw_nfft);
  h5_write(archive, "M4_iw", S.M4_iw);
  h5_write(archive, "M4pp_iw", S.M4pp_iw);
  h5_write(archive, "M4ph_iw", S.M4ph_iw);
  h5_write(archive, "M3pp_iw_nfft", S.M3pp_iw_nfft);
  h5_write(archive, "M3ph_iw_nfft", S.M3ph_iw_nfft);
  h5_write(archive, "M3pp_tau", S.M3pp_tau);
  h5_write(archive, "M3ph_tau", S.M3ph_tau);
  h5_write(archive, "M3xph_tau", S.M3xph_tau);
  h5_write(archive, "chi2pp_tau", S.chi2pp_tau);
  h5_write(archive, "chi2ph_iw", S.chi2ph_tau);
  h5_write(archive, "chiAB_tau", S.chiAB_tau);
}

MAKE_MAIN;
