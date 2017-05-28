from itertools import product
import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from pytriqs.utility.h5diff import h5diff
from triqs_ctint import SolverCore

test_name = "densdens"

######## physical parameters ########
U=0.5
mu=U/4.
beta = 10.0

######## simulation parameters ########
n_cyc = 1000

# --------- set up static interactions and the block structure ---------
N_states=1
block_names = ['up','dn']
composite_blocks = [ b[0]+'|'+b[1] for b in product(block_names,block_names)]
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
               n_tau = 100001,
               use_D = True,
               n_tau_dynamical_interactions = 1025,
               n_iw_dynamical_interactions = 200)

# --------- Initialize the non-interacting Green's function ----------
semicirc = S.G0_iw[block_names[0]].copy()
semicirc << SemiCircular(1.0)
for n,g in S.G0_iw: g << inverse(iOmega_n + mu - 1.0 * semicirc)

# --------- Set up time-dependent interactions ----------
# Boson Frequency
w0=1.0

# Dynamic Density-Density Interaction
D = 0.5;
for b1 in block_names:
  for b2 in block_names:
    for i in range(N_states):
      for j in range(N_states):
        if i==j and b1==b2:
          S.D0_iw(b1,b2)[i,j]  << D**2*(inverse(iOmega_n-w0)-inverse(iOmega_n+w0))  
        else:
          S.D0_iw(b1,b2)[i,j]  << D**2*(inverse(iOmega_n-w0)-inverse(iOmega_n+w0))
 
# --------- Solve! ----------
S.solve(h_int=h_int,
        alpha = alpha,
        n_cycles = n_cyc,
        length_cycle = 50,
        n_warmup_cycles = 0,
        random_seed = 34788,
        measure_M_tau = True,
        measure_M_iw = True,
        measure_F_tau = True,
        post_process = False )

# -------- Save in archive ---------
A = HDFArchive("%s.out.h5"%test_name,'w')
A["G0_iw"] = S.G0_iw
A["M_tau"] = S.M_tau
A["M_iw_nfft"] = S.M_iw_nfft
A["F_tau"] = S.F_tau

# -------- Compare ---------
h5diff("%s.out.h5"%test_name, "%s.ref.h5"%test_name)
