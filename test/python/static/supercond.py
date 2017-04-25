from itertools import product
from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from numpy import array, zeros, identity, matrix


test_name = 'supercond'

#superconducting calc: U n_up n_down - mu n = Utilde n_1 n_2 - mutilde (n1+n2)

# set up a few parameters
U = -1.0
Utilde=-U
beta = 20.0
mu = matrix([[1.3-U,0],[0,-1.3]])
eps = matrix([[0.2,0.1],[0.1,-0.2]])
delta = 0.35

block_names = ['up']
composite_blocks = [ b[0]+'|'+b[1] for b in product(block_names,block_names)]

h_int = Utilde * n(block_names[0],0)*n(block_names[0],1)

# alpha[block][index,s]
diag = 0.5 + delta
odiag = 0.5 - delta
alpha = [ [[diag,odiag],[diag, odiag]] ]

gf_struct = {block_names[0] : [0,1]}

# Construct the segment solver
from ctint import SolverCore
S = SolverCore(beta = beta, 
               gf_struct = gf_struct,
               n_iw = 200,  
               n_tau = 100001 )

# Initialize the Green's function
for n,g0 in S.G0_iw:
  g0 << inverse(iOmega_n + mu - inverse(iOmega_n - eps));

# Solve!
S.solve(h_int=h_int,
        alpha = alpha,
        n_cycles = 1000,
        length_cycle = 50,
        n_warmup_cycles = 0,
        random_seed = 34788,
        measure_M_tau = True,
        measure_M_iw = True,
        measure_F_tau = True,
        post_process = False )

# Save in archive
A = HDFArchive("%s.out.h5"%test_name,'w')
A["G0_iw"] = S.G0_iw
A["M_tau"] = S.M_tau
A["M_iw_nfft"] = S.M_iw_nfft
A["F_tau"] = S.F_tau

from pytriqs.utility.h5diff import h5diff
h5diff("%s.out.h5"%test_name, "%s.ref.h5"%test_name)
