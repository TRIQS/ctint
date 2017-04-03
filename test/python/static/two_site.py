from itertools import product
from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.archive import *
from numpy import array, zeros, identity, matrix

test_name = 'two_site'

# set up a few parameters
U = 1.0
V= 0.3
beta = 20.0
mu = 1.3
eps = matrix([[0.2,0.1],[0.1,0.2]])
delta = 0.35

block_names = ['up','down']
composite_blocks = [ b[0]+'|'+b[1] for b in product(block_names,block_names)]

h_int =  U * (   n('up',0)*n('down',0)\
               + n('up',1)*n('down',1) )\
       + V * (   n('up',0)*n('up',1)\
               + n('down',0)*n('down',1)\
               + n('up',0)*n('down',1)\
               + n('up',1)*n('down',0) )

# alpha[block][index,s]
diag = 0.5 + delta
odiag = 0.5 - delta
alpha = [ [[diag,odiag],[diag, odiag]], [[odiag,diag], [odiag, diag] ]]

gf_struct = dict.fromkeys(block_names, [0,1])

# Construct the segment solver
from ctint import SolverCore
S = SolverCore(beta = beta, 
               gf_struct = gf_struct,
               n_iw = 200,  
               n_tau = 100001)

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
        measure_F_tau = True,
        post_process = False )

# Save in archive
A = HDFArchive("%s.out.h5"%test_name,'w')
A["G0_iw"] = S.G0_iw
A["M_tau"] = S.M_tau
A["F_tau"] = S.F_tau

from pytriqs.utility.h5diff import h5diff
h5diff("%s.out.h5"%test_name, "%s.ref.h5"%test_name)
