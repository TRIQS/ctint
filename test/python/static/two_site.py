from triqs_ctint import Solver

from itertools import product
import triqs.utility.mpi as mpi
from triqs.gf import *
from h5 import *
from triqs.operators import *
from triqs.utility.h5diff import h5diff
from triqs.operators.util.hamiltonians import h_int_kanamori
from numpy import matrix, array

test_name = 'two_site'

######## physical parameters ########
U = 1.                          # Density-density interaction for opposite spins
Up = 1.                         # Density-density interaction for equal spins
J = 0.2                         # Hunds coupling
mu = U/4.0
beta = 10.0
eps = matrix([  [0.2,0.1],
                [0.1,0.2]])

######## simulation parameters ########
n_cyc = 1000

# --------- set up static interactions and the block structure ---------
block_names = ['dn','up']
orb_names = [0,1]
gf_struct = [(bl, len(orb_names)) for bl in block_names]
h_int = h_int_kanamori(block_names, orb_names,
                        array([[0,      Up-3*J ],
                               [Up-3*J, 0      ]]), # Interaction for equal spins
                        array([[U,      U-2*J  ],
                               [U-2*J,  U      ]]),   # Interaction for opposite spins
                        J,off_diag=True)

# --------- Construct the ctint solver ----------
S = Solver(beta = beta,
               gf_struct = gf_struct,
               n_iw = 200,
               n_tau = 100001)

# --------- Initialize the non-interacting Green's function ----------
for bl, g_bl in S.G0_iw: g_bl << inverse(iOmega_n + mu - inverse(iOmega_n - eps));

# --------- Solve! ----------
S.solve(h_int=h_int,
        n_cycles = n_cyc,
        length_cycle = 100,
        n_warmup_cycles = 100,
        random_seed = 34788,
        measure_histogram = True,
        measure_density = True,
        measure_M4_iw = True,
        n_iw_M4 = 5,
        nfft_buf_size = 50,
        measure_M3pp_tau = True,
        measure_M3ph_tau = True,
        measure_M3xph_tau = True,
        n_iw_M3 = 10,
        n_iW_M3 = 10,
        n_tau_M3 = 41,
        measure_chi2pp_tau = True,
        measure_chi2ph_tau = True,
        n_iw_chi2 = 10,
        n_tau_chi2 = 21,
        measure_chiAB_tau = True,
        chi_A_vec = [n('up',0) + n('dn', 0)],
        chi_B_vec = [n('up',0) + n('dn', 0)],
        post_process = True )

# -------- Save in archive ---------
with HDFArchive("%s.out.h5"%test_name,'w') as arch:
    arch["G0_iw"] = S.G0_iw
    arch["G_iw"] = S.G_iw
    arch["G2_iw"] = S.G2_iw
    arch["chi3pp_iw"] = S.chi3pp_iw
    arch["chi3ph_iw"] = S.chi3ph_iw
    arch["chi3xph_iw"] = S.chi3xph_iw
    arch["chi2pp_iw"] = S.chi2pp_iw
    arch["chi2ph_iw"] = S.chi2ph_iw
    arch["chiAB_iw"] = S.chiAB_iw
    arch["chi2pp_tau_from_M3"] = S.chi2pp_tau_from_M3
    arch["chi2ph_tau_from_M3"] = S.chi2ph_tau_from_M3
    arch["chi2xph_tau_from_M3"] = S.chi2xph_tau_from_M3

# -------- Compare ---------
h5diff("%s.out.h5"%test_name, "%s.ref.h5"%test_name)
