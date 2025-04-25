#include <triqs_ctint/solver_core.hpp>

#include <h5/h5.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

using namespace triqs_ctint;

int main() { // NOLINT
  // --------- physical parameters ----------
  double const U    = 1.0;     // Density-density interaction
  double const t    = 1.0;     // Hopping
  double const mu   = U / 2.0; // Chemical Potential
  double const beta = 100;     // Inverse temperature

  // --------- simulation parameters ----------
  int const n_cyc = 500;

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

  // --------- Solve! ----------
  long n_terms = std::distance(h_int.begin(), h_int.end());
  alpha_t alpha(n_terms, 2, 2, 1);
  double const delta = 0.1;
  for (long l = 0; l < n_terms; ++l) { alpha(l, range::all, range::all, 0) = nda::matrix<double>{{0.5 - delta, 0.}, {0., 0.5 + delta}}; };

  solve_params_t ps;
  ps.h_int           = h_int;
  ps.alpha           = alpha;
  ps.n_cycles        = n_cyc;
  ps.length_cycle    = 50;
  ps.n_warmup_cycles = 100;
  ps.random_seed     = 34788;
  ps.nfft_buf_size   = 100000;
  ps.measure_M3pp_iw = true;
  ps.measure_M3ph_iw = true;
  ps.n_iw_M3         = 24;
  ps.n_iW_M3         = 24;
  ps.post_process    = false;

  S.solve(ps);
}
