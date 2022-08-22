from numpy import zeros,matrix, array,sinh,cosh, cos, sin, exp, arctan, linspace
import triqs.utility.mpi as mpi
from triqs.gf import *
from h5 import *
from triqs.operators import *
from cthyb_spin import SolverCore

n_cycles =10000000 # 40000
max_time=1800
n_iw = 32 #number of freqs in vertex

U=1.0
mu=U/2.
beta = 20.0
block_names = ['up','down']
h_int = U * n(block_names[0],0)*n(block_names[1],0)
gf_struct = {block_names[0] : [0], block_names[1] : [0]}

S = SolverCore(beta = beta, gf_struct = gf_struct,
                 n_w = 200,
                 n_w_b_nn = 200,
                 n_tau_nn = 401,
                 )
semicirc = S.G0_iw[block_names[0]].copy()
semicirc << SemiCircular(1.0)
for n,g in S.G0_iw: g << inverse(iOmega_n + mu - 1.0 * semicirc)

S.solve(h_int=h_int,
        max_time = max_time,
        n_cycles = n_cycles,
        length_cycle = 200,
        n_warmup_cycles = 5000,
        measure_gw = True,        measure_g2w = True,
        measure_g3w = True,
        measure_ft = True,        measure_nnt = True,
        measure_gt = True,
        measure_nnw = True,
        measure_nn = True,
        measure_hist = True,
        n_w_f_vertex = n_iw,
        n_w_b_vertex = n_iw,
        )

if mpi.is_master_node():
 A = HDFArchive("cthyb.out.h5",'w')
 A["G_iw"] = S.G_iw
 A["G_tau"] = S.G_tau
 A["F_tau"] = S.F_tau
 A["nn_tau"] = S.nn_tau
 A["nn_iw"] = S.nn_iw
 A["nn"] = S.nn
 A["G_2w"] = S.G_2w
 A["G_3w"] = S.G_3w

