#include "solver_core.hpp"
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs_ctint;

TEST(CtInt, Anderson) {
  triqs::mpi::communicator world;

  // System Parameters
  double mu      = 1.3;
  double epsilon = 0.2;
  double delta   = 0.35;
  double U       = 1.0;

  // Construct Parameters
  constr_params_t pc;
  pc.n_tau     = 10000;
  pc.n_iw      = 500;
  pc.beta      = 20.0;
  pc.gf_struct = {{"up", {0}}, {"down", {0}}};
  pc.use_D     = false;
  pc.use_Jperp = false;

  // Solve Parameters
  solve_params_t ps;
  ps.h_int     = U * n("up", 0) * n("down", 0);
  ps.use_alpha = true;
  ps.n_s       = 2;
  // alpha parameters
  typedef std::vector<array<double, 2>> alpha_t;
  alpha_t alpha(2);
  double diag             = 0.5 + delta;
  double odiag            = 0.5 - delta;
  alpha[0]                = array<double, 2>{{diag, odiag}};
  alpha[1]                = array<double, 2>{{odiag, diag}};
  ps.alpha                = alpha;
  ps.n_cycles             = 1000;
  ps.length_cycle         = 50;
  ps.n_warmup_cycles      = 1000;
  ps.random_seed          = 34788;
  ps.random_name          = "";
  ps.use_double_insertion = false;
  ps.max_time             = -1;
  ps.verbosity            = 3;
  ps.measure_average_sign = true;
  ps.measure_M_tau        = true;
  ps.measure_M_iw         = false;
  ps.nfft_buf_size        = 10000;
  ps.post_process         = false;

  solver_core ctqmc(pc);

  auto sha = triqs::arrays::make_shape(1, 1);

  // set delta(i\omega_n) and g0(i\omega_n)
  auto delta_w = gf<imfreq, matrix_valued>{{pc.beta, Fermion, pc.n_iw}, sha};
  delta_w(iw_) << 1.0 / (iw_ - epsilon);
  ctqmc.G0_iw()[0](iw_) << 1.0 / (iw_ + mu + (-1.0) * delta_w(iw_));
  ctqmc.G0_iw()[1](iw_) << 1.0 / (iw_ + mu + (-1.0) * delta_w(iw_));

  // Solve!
  ctqmc.solve(ps);

  // Save the results
  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.out.h5", 'w');
    h5_write(G_file, "M_tau", ctqmc.M_tau);
  }

  // Compare
  if (world.rank() == 0) {
    triqs::h5::file G_file("anderson.ref.h5", 'r');
    block_gf<imtime, matrix_valued> M_tau_ref;
    h5_read(G_file, "M_tau", M_tau_ref);
    EXPECT_BLOCK_GF_NEAR(M_tau_ref, *ctqmc.M_tau);
  }
}
MAKE_MAIN;
