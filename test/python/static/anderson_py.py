from itertools import product
import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from pytriqs.utility.h5diff import h5diff
from triqs_ctint import SolverCore

test_name = "anderson_py"

######## physical parameters ########
U = 1.0
mu = U/4.0
beta = 10.0
eps = 0.2

######## simulation parameters ########
n_cyc = 1000

# --------- set up static interactions and the block structure ---------
block_names = ['up','dn']
gf_struct = dict.fromkeys(block_names, [0])
h_int = U * n(block_names[0],0)*n(block_names[1],0)

# --------- Define alpha tensor ---------
delta = 0.1
diag = 0.5 + delta
odiag = 0.5 - delta
alpha = [ [[diag,odiag]], [[odiag,diag]] ] # alpha[block][index,s]

# --------- Construct the ctint solver ----------
S = SolverCore(beta = beta,
               gf_struct = gf_struct,
               n_iw = 200,
               n_tau = 100001)

# --------- Initialize the non-interacting Green's function ----------
for n,g0 in S.G0_iw:
  g0 << inverse(iOmega_n + mu - inverse(iOmega_n - eps));

# --------- Solve! ----------
S.solve(h_int=h_int,
        alpha = alpha,
        n_cycles = n_cyc,
        length_cycle = 50,
        n_warmup_cycles = 0,
        random_seed = 34788,
        measure_M_tau = True,
        measure_M4_iw = True,
        n_iw_M4 = 5,
        nfft_buf_size = 50,
        measure_M3pp_tau = True,
        measure_M3ph_tau = True,
        n_iw_M3 = 10,
        n_iW_M3 = 10,
        n_tau_M3 = 100,
        measure_chi2pp_tau = True,
        measure_chi2ph_tau = True,
        n_iw_chi2 = 10,
        n_tau_chi2 = 30,
        post_process = True )

# -------- Save in archive ---------
with HDFArchive("%s.out.h5"%test_name,'w') as arch:
    arch["G0_iw"] = S.G0_iw
    arch["G_iw"] = S.G_iw
    arch["G2_iw"] = S.G2_iw
    arch["chi3pp_iw"] = S.chi3pp_iw
    arch["chi3ph_iw"] = S.chi3ph_iw
    arch["chi2pp_iw"] = S.chi2pp_iw
    arch["chi2ph_iw"] = S.chi2ph_iw

# -------- Save Solver and Retrieve and solve
with HDFArchive("solver.h5",'w') as arch:
    arch["S"] = S

with HDFArchive("solver.h5",'r') as arch:
    S = arch["S"]
    S.solve(**S.solve_params)

# -------- Save 2nd run in archive ---------
with HDFArchive("%s.out_2nd.h5"%test_name,'w') as arch:
    arch["G0_iw"] = S.G0_iw
    arch["G_iw"] = S.G_iw
    arch["G2_iw"] = S.G2_iw
    arch["chi3pp_iw"] = S.chi3pp_iw
    arch["chi3ph_iw"] = S.chi3ph_iw
    arch["chi2pp_iw"] = S.chi2pp_iw
    arch["chi2ph_iw"] = S.chi2ph_iw

# -------- Compare ---------
h5diff("%s.out.h5"%test_name, "%s.ref.h5"%test_name)
h5diff("%s.out.h5"%test_name, "%s.out_2nd.h5"%test_name)
