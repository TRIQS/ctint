#include "solver_core.hpp"
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs_ctint;

TEST(CtInt, Anderson) {
  triqs::mpi::communicator world;

  // System Parameters
  double mu    = 0.0;
  double delta = 0.35;
  double U     = 1.0;

  // Discrete bath energies and hoppings
  std::vector<double> energ = {-0.7, -0.15, 0.15, 0.7}; // Vertex Paper ( optimized for Beta=20 ) to fit box DOS
  std::vector<double> hopp  = {0.45, 0.34, 0.34, 0.45};

  // The alpha function
  typedef std::vector<array<double, 2>> alpha_t;
  alpha_t alpha(2);
  double diag  = 0.5 + delta;
  double odiag = 0.5 - delta;
  alpha[0]     = array<double, 2>{{diag, odiag}};
  alpha[1]     = array<double, 2>{{odiag, diag}};

  // Construct Parameters
  constr_params_t pc;
  pc.beta      = 10.0;
  pc.gf_struct = {{"up", {0}}, {"down", {0}}};
  pc.n_tau     = 10000;
  pc.n_iw      = 500;
  pc.use_D     = false;
  pc.use_Jperp = false;

  // Solve Parameters
  solve_params_t ps;
  ps.h_int                = U * n("up", 0) * n("down", 0);
  ps.use_alpha            = true;
  ps.n_s                  = 2;
  ps.alpha                = alpha;
  ps.n_cycles             = 1e+6;
  ps.length_cycle         = 50;
  ps.n_warmup_cycles      = 5000;
  ps.use_double_insertion = false;
  ps.max_time             = -1;
  ps.verbosity            = 3;
  ps.post_process         = true;
  //ps.random_seed          = 34788 + mpi..;
  //ps.random_name          = "";

  ps.measure_average_sign = true;

  ps.measure_M_tau = false;
  ps.measure_M_iw  = false;

  //ps.measure_M4_iw  = true;
  //ps.n_iw_M4        = 32;

  //ps.measure_M3pp_iw  = true;
  //ps.measure_M3ph_iw  = true;
  //ps.measure_M3xph_iw = true;
  //ps.n_iw_M3          = 64;

  //ps.measure_M2pp_tau  = true;
  //ps.measure_M2ph_tau  = true;
  //ps.measure_M2xph_tau = true;
  //ps.n_tau_M2          = 1000;
  //ps.n_iw_M2           = 128;

  ps.nfft_buf_size = 100000;

  solver_core ctqmc(pc);

  // set up hybridization delta(i\omega_n)
  auto delta_w = gf<imfreq, matrix_valued>{{pc.beta, Fermion, pc.n_iw}, make_shape(1, 1)};
  delta_w()    = 0;
  for (int i : range(energ.size())) {
    auto term = gf<imfreq, matrix_valued>{{pc.beta, Fermion, pc.n_iw}, make_shape(1, 1)};
    term(iw_) << hopp[i] * hopp[i] / (iw_ - energ[i]);
    delta_w += term;
  }

  // set up non-interacting Green function g0(i\omega_n)
  ctqmc.G0_iw()[0](iw_) << 1.0 / (iw_ + mu + (-1.0) * delta_w(iw_));
  ctqmc.G0_iw()[1](iw_) << 1.0 / (iw_ + mu + (-1.0) * delta_w(iw_));

  // Solve!
  ctqmc.solve(ps);
}

MAKE_MAIN;
