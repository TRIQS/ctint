from triqs_ctint import Solver

from triqs.gf import *
from h5 import HDFArchive
from triqs.utility import mpi
from triqs.utility.comparison_tests import *
from triqs.operators import *

U    = 1.0
mu   = U/2.0
beta = 10

# Create Solver
cp = {
  "beta": beta,
  "gf_struct" : [["up", 1], ["down", 1]]
}
S = Solver(**cp)
for bl, g_bl in S.G0_iw: g_bl << inverse(iOmega_n + mu);

# Parameters for the Run
sp = {
  "h_int"           : U * n("up", 0) * n("down", 0),
  "n_s"             : 1,
  "alpha"           : [[[0.5 + 0.1]], [[0.5 - 0.1]]],
  "length_cycle"    : 50,
  "n_warmup_cycles" : 1000,
  "n_cycles"        : 1000
}

# Solve the model
S.solve(**sp)

if mpi.is_master_node(): HDFArchive('Solver.h5', 'w')['S'] = S

# Rerun the Solver and Compare
S_old = HDFArchive('Solver.h5', 'r')['S']
S_old.solve(**S_old.last_solve_params)
assert_block_gfs_are_close(S_old.G_iw, S.G_iw, 1e-14)
