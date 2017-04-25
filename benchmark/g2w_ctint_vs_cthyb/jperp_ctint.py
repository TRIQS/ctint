from numpy import zeros,matrix, array,sinh,cosh, cos, sin, exp, arctan, linspace
import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
from ctint import SolverCore

n_cycles =10000000 # 40000
max_time=1800
n_iw = 64 #number of frequencies in vertex

n_tau_M4t=200

U=1.0
mu=U/2.
beta = 20.0
block_names = ['up','down']
gf_struct = {block_names[0] : [0], block_names[1] : [0]}
h_int = U * n(block_names[0],0)*n(block_names[1],0)
delta = 0.3
diag = 0.5 + delta
odiag = 0.5 - delta
alpha = [ [[diag,odiag]], [[odiag,diag]] ] # alpha[block][index,s]

S = SolverCore(beta = beta, gf_struct = gf_struct, n_iw = 200,  n_tau_g0  = 1025,
                 n_tau_f  = 1025,
                 n_tau_dynamical_interactions = 1025,
                 n_iw_dynamical_interactions = 200,
                 n_tau_M4t = n_tau_M4t,
                 n_w_f_M4w = n_iw, 
                 n_w_b_M4w = n_iw,
                 n_tau_g2t = 2*n_iw+1,
                 n_w_f_g2w = n_iw,
                 n_w_b_g2w = n_iw,
                 n_tau_nnt = 401,
                 )
semicirc = S.G0_iw[block_names[0]].copy()
semicirc << SemiCircular(1.0)
for n,g in S.G0_iw: g << inverse(iOmega_n + mu - 1.0 * semicirc)

S.solve(h_int=h_int,
        alpha = alpha,
        max_time = max_time,
        n_cycles = n_cycles,  
        length_cycle = 200, 
        n_warmup_cycles = 5000,
        only_sign = False,
        measure_gw = True,        measure_Mt = True,
        measure_ft = True,        measure_nnt = True,
        measure_hist = True,   
        measure_nn = True,
        measure_M4t = False,
        measure_g2t = True,
        print_post_process = False,     g2t_indep = [],
        )

if mpi.is_master_node():
 A = HDFArchive("ctint_n_tau_M4t_%s_niw_%s_g2t.out.h5"%(n_tau_M4t, n_iw),'w')
 A["G_iw"] = S.G_iw
 A["M_tau"] = S.M_tau
 A["F_tau"] = S.F_tau
 A["nn_tau"] = S.nn_tau
 A["nn_iw"] = S.nn_iw
 A["nn"] = S.nn
 A["G_2w"] = S.chi3_iwinu
 A["hist"] = S.histogram

