#include <triqs_ctint/solver_core.hpp>

#include <triqs/h5.hpp>
#include <mpi/mpi.hpp>

#include <string>

using namespace triqs_ctint;

std::string file_name = "Solver_old.h5";

int main(int argc, char **argv) {

  auto env  = mpi::environment(argc, argv);
  auto comm = mpi::communicator();

  // File containing a stored triqs_ctint solver object

  auto arch  = triqs::h5::file(file_name, 'r');
  auto S_old = solver_core::h5_read_construct(arch, "S");
  auto cp    = S_old.constr_params;
  auto sp    = S_old.last_solve_params.value();

  auto S    = solver_core{cp};
  S.G0_iw() = S_old.G0_iw;

  S.solve(sp);

  if (comm.rank() == 0) {
    auto arch_out = triqs::h5::file("Solver.h5", 'w');
    h5_write(arch_out, "S", S);
  }
}
