from itertools import product
import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from pytriqs.utility.h5diff import h5diff
from triqs_ctint import SolverCore
from numpy import matrix

test_name = 'two_site'

######## physical parameters ########
U = 1.0
V = 0.3
mu = U/4.0
beta = 10.0
eps = matrix([[0.2,0.1],[0.1,0.2]])

######## simulation parameters ########
n_cyc = 1000

# --------- set up static interactions and the block structure ---------
block_names = ['up','dn']
gf_struct = dict.fromkeys(block_names, [0,1])
h_int =  U * ( n('up',0)*n('dn',0) + n('up',1)*n('dn',1) ) \
       + V * ( n('up',0)*n('up',1) + n('dn',0)*n('dn',1) + \
               n('up',0)*n('dn',1) + n('dn',0)*n('up',1) )

# --------- Define alpha tensor ---------
delta = 0.1
diag = 0.5 + delta
odiag = 0.5 - delta
alpha = [ [[diag,odiag],[diag,odiag]], [[odiag,diag],[odiag,diag]] ] # alpha[block][index,s]

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
        measure_M_iw = True,
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
A = HDFArchive("%s.out.h5"%test_name,'w')
A["G0_iw"] = S.G0_iw
A["G_iw"] = S.G_iw
A["M_iw_nfft"] = S.M_iw_nfft
A["chi3pp_iw"] = S.chi3pp_iw
A["chi3ph_iw"] = S.chi3ph_iw
A["chi2pp_iw"] = S.chi2pp_iw
A["chi2ph_iw"] = S.chi2ph_iw

# -------- Compare ---------
h5diff("%s.out.h5"%test_name, "%s.ref.h5"%test_name)
