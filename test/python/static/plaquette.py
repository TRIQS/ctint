from triqs_ctint import Solver

from itertools import product
import triqs.utility.mpi as mpi
from triqs.gfs import *
from h5 import *
from triqs.operators import c, c_dag, n
from triqs.utility.h5diff import h5diff
from numpy import matrix, array, zeros, reshape, eye

test_name = "plaquette"

######## physical parameters ########
U = 1.0  # Density-density interaction
t = 1.0  # Hopping
mu = U / 2.0  # Chemical Potential
beta = 2.0  # Inverse temperature

######## simulation parameters ########
n_cyc = 50
# NFFT transform accuracy, set far below the fixed h5diff tolerance used at the end.
# The machine-dependent NFFT auto-dispatch can shift NFFT-derived
# observables (M4_iw, M3*, chi2*, chi3*/G2) by a value-dependent multiple of nfft_tol between runs; keeping
# nfft_tol well under the comparison tolerance keeps that drift far below it.
nfft_tol = 1e-14

# --------- Define hopping matrix and interaction hamiltonian ----------

Nx = 2
Ny = 2
n_orb = Nx * Ny

hloc0_mat = zeros((n_orb, n_orb))

# diagonal terms
hloc0_mat -= eye(n_orb)*mu

hloc0_arr = reshape(hloc0_mat, (Nx, Ny, Nx, Ny))
for x in range(Nx - 1):
    hloc0_arr[x, :, x + 1, :] -= t*eye(Ny)
    hloc0_arr[x + 1, :, x, :] -= t*eye(Ny)
for y in range(Ny - 1):
    hloc0_arr[:, y, :, y + 1] -= t*eye(Nx)
    hloc0_arr[:, y + 1, :, y] -= t*eye(Nx)

h_int = sum(U * n("up", i) * n("dn", i) for i in range(n_orb))

# --------- set up block structure ---------

block_names = ["dn", "up"]
gf_struct = [(bl, n_orb) for bl in block_names]

# --------- Construct the ctint solver ----------
S = Solver(beta = beta,
               gf_struct = gf_struct,
               n_tau = 201,
               dlr_wmax = 10.0,
               dlr_eps = 1e-7)

# --------- Initialize the non-interacting Green's function ----------
for bl, g_bl in S.G0_iw: g_bl << inverse(iOmega_n - hloc0_mat)

# --------- Solve! ----------
S.solve(h_int=h_int,
        n_cycles = n_cyc,
        length_cycle = 100,
        n_warmup_cycles = 100,
        random_seed = 34788,
        measure_histogram = True,
        measure_density = True,
        measure_M_tau = True,
        measure_M_iw = True,
        measure_M4_iw = True,
        measure_M4pp_iw = True,
        measure_M4ph_iw = True,
        n_iw_M4 = 2,
        n_iW_M4 = 2,
        nfft_buf_size = 100000,
        nfft_tol = nfft_tol,
        measure_M3pp_iw = True,
        measure_M3ph_iw = True,
        measure_M3pp_iw_full = True,
        measure_M3ph_iw_full = True,
        measure_M3pp_tau = True,
        measure_M3ph_tau = True,
        measure_M3xph_tau = True,
        n_iw_M3 = 4,
        n_iW_M3 = 4,
        n_tau_M3 = 4,
        measure_chi2pp_tau = True,
        measure_chi2ph_tau = True,
        measure_chiAB_tau = True,
        chi_A_vec = [n('up',0) + n('dn', 0)],
        chi_B_vec = [n('up',0) + n('dn', 0)],
        post_process = False,
        use_double_insertion = True)

# -------- Save in archive ---------
with HDFArchive("%s.out.h5"%test_name,'w') as arch:
    arch["histogram"] = S.histogram
    arch["density"] = S.density
    arch["M_tau"] = S.M_tau
    arch["M_iw"] = S.M_iw_nfft
    arch["M4_iw"] = S.M4_iw
    arch["M4pp_iw"] = S.M4pp_iw
    arch["M4ph_iw"] = S.M4ph_iw
    # M3pp_iw_nfft and M3ph_iw_nfft live on DLR2D meshes whose construction
    # is currently BLAS-implementation dependent. This may be changed in the
    # future. For now we exclude these quantities from the test comparison.
    # The tau variants already cover these measurements.
    #arch["M3pp_iw_nfft"] = S.M3pp_iw_nfft
    #arch["M3ph_iw_nfft"] = S.M3ph_iw_nfft
    #arch["M3pp_iw_nfft_full"] = S.M3pp_iw_nfft_full
    #arch["M3ph_iw_nfft_full"] = S.M3ph_iw_nfft_full
    arch["M3pp_tau"] = S.M3pp_tau
    arch["M3ph_tau"] = S.M3ph_tau
    arch["M3xph_tau"] = S.M3xph_tau
    arch["chi2pp_tau"] = S.chi2pp_tau
    arch["chi2ph_iw"] = S.chi2ph_tau
    arch["chiAB_tau"] = S.chiAB_tau

# -------- Compare ---------
h5diff("%s.out.h5"%test_name, "%s.ref.h5"%test_name, precision=2e-5)
