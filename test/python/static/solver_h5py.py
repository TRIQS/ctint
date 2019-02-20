from triqs_ctint import Solver

from pytriqs.gf import *
from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pytriqs.utility.comparison_tests import *
from pytriqs.operators import *

U    = 1.0
mu   = U/2.0
beta = 10

# Create Solver
cp = {
  "beta": beta,
  "gf_struct" : [["up", [0]], ["down", [0]]]
}
S = Solver(**cp)
for bl, g_bl in S.G0_iw: g_bl << inverse(iOmega_n + mu);

# Parameters for the Run
sp = {
  "h_int"              : U * n("up", 0) * n("down", 0),
  "length_cycle"       : 50,
  "n_warmup_cycles"    : 1000,
  "n_cycles"           : 1000,
  "measure_M3ph_tau"   : True,
  "n_iw_M3"            : 10,
  "n_iW_M3"            : 10,
  "n_tau_M3"           : 121,
  "measure_chi2ph_tau" : True,
  "n_iw_chi2"          : 10,
  "n_tau_chi2"         : 100
}

# Solve the model
S.solve(**sp)

if mpi.is_master_node(): HDFArchive('Solver.h5', 'w')['S'] = S

# Rerun the Solver and Compare
S2 = HDFArchive('Solver.h5', 'r')['S']
S2.solve(**S2.last_solve_params)

assert_block_gfs_are_close(S.G_iw, S2.G_iw, 1e-16)
assert_block2_gfs_are_close(S.M3ph_tau, S2.M3ph_tau, 1e-16)
assert_block2_gfs_are_close(S.chi2ph_tau, S2.chi2ph_tau, 1e-16)
