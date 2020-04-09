#include <triqs_ctint/solver_core.hpp>

#include <h5/h5.hpp>
#include <mpi/mpi.hpp>
#include <triqs/gfs/gf_tests.hpp>
#include <triqs/clef.hpp>

#include <string>
#include <vector>

using namespace triqs_ctint;
using namespace h5;

int main(int argc, char **argv) {

  auto env  = mpi::environment(argc, argv);
  auto comm = mpi::communicator();

  // Hubbard Atom parameters
  double U    = 1.0;
  double mu   = U/2.0;
  double beta = 10;

  // Create Solver
  constr_params_t cp{};
  cp.beta      = beta;
  cp.gf_struct = {{"up", {0}}, {"down", {0}}};
  auto S       = triqs_ctint::solver_core{cp};
  S.G0_iw[0][iw_] << 1.0 / (iw_ + mu);
  S.G0_iw[1][iw_] << 1.0 / (iw_ + mu);

  // Parameters for the Run
  solve_params_t sp;
  sp.h_int           = U * n("up", 0) * n("down", 0);
  sp.alpha           = {{{0.5 + 0.1}}, {{0.5 - 0.1}}};
  sp.length_cycle    = 50;
  sp.n_warmup_cycles = 1000;
  sp.n_cycles        = 1000;

  // Solve the model
  S.solve(sp);

  // Store the Solver
  if (comm.rank() == 0) h5_write(file("Solver.h5", 'w'), "S", S);

  // Rerun the Solver and Compare
  auto S_old = solver_core::h5_read_construct(file("Solver.h5", 'r'), "S");
  S_old.solve(S_old.last_solve_params.value());
  assert_block_gfs_are_close(S_old.G_iw, S.G_iw, 1e-15);
}
