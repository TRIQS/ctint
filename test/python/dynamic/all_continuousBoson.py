from itertools import product
from pytriqs.gf.local import *
from pytriqs.archive import *
from pytriqs.operators import *
from pytriqs.utility.h5diff import h5diff
from pytriqs.gf.local.descriptors import Function
from numpy import linspace
from ctint import SolverCore

test_name = 'all_continuousBoson'

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
               use_Jperp = True,
               n_tau_dynamical_interactions = 1025,
               n_iw_dynamical_interactions = 200)

# --------- Initialize the non-interacting Green's function ----------
semicirc = S.G0_iw[block_names[0]].copy()
semicirc << SemiCircular(1.0)
for n,g in S.G0_iw: g << inverse(iOmega_n + mu - 1.0 * semicirc)

# --------- Set up time-dependent interactions ----------
w0=1.0
s=0.2
nu = lambda n,beta : (2*n*pi/beta) #bosonic frequences
Eps=linspace(0,w0,2500) #eps grid
Dz_expr = lambda w, g : (-1.)*(.5*g**2) *((s+1)*w0**(-s-1))*sum([eps**(s+1)/(w.imag*w.imag+eps**2) for eps in Eps])*(Eps[1]-Eps[0]) #Dz away from 0
Dz0 = lambda g : (-1.)*(.5*g**2)*(s+1.)/s/w0 #Dz at zero (Dz_expr ill-defined at w=0)
Dz_expr_reg = lambda w, g: Dz_expr(w,g) if abs(w)>0.0 else Dz0(g) #Dz for any w

# Dynamic Spin-Spin Interaction
J = 0.5;
for i in range(N_states):
  for j in range(N_states):
    S.Jperp_iw[i,j] << Function(lambda w : 4.0*Dz_expr_reg(w, J)) 
   
# Dynamic Density-Density Interaction
D = 0.5;
for b1 in block_names:
  for b2 in block_names:
    for i in range(N_states):
      for j in range(N_states):
        if i==j and b1==b2:
          S.D0_iw[b1+"|"+b2][i,j] << Function(lambda w : Dz_expr_reg(w, D))
        else:
          S.D0_iw[b1+"|"+b2][i,j] << Function(lambda w : Dz_expr_reg(w, D))

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
