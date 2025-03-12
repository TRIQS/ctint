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

  auto const hloc0_mat = -nda::matrix<double>{
     {mu, t, 0, t},
     {t, mu, t, 0},
     {0, t, mu, t},
     {t, 0, t, mu},
  };

  many_body_operator h_int = U * n("up", 0) * n("dn", 0) + U * n("up", 1) * n("dn", 1) + U * n("up", 2) * n("dn", 2) + U * n("up", 3) * n("dn", 3);

  // --------- set up block structure ---------

  gf_struct_t gf_struct{{"dn", 4}, {"up", 4}};

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

  // --------- Solve! ----------
  double const delta  = 0.1;
  double const a_upup = 0.5 - delta;
  double const a_dndn = 0.5 + delta;
  alpha_t alpha(4, 2, 2, 1);
  for (long i = 0; i < 4; ++i) { alpha(i, range::all, range::all, 0) = nda::matrix<double>{{a_upup, 0.}, {0., a_dndn}}; };

  solve_params_t ps;
  ps.alpha              = alpha;
  ps.h_int              = h_int;
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
  ps.nfft_buf_size      = 50;
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
