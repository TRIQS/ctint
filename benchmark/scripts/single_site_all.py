from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.operators import *
import numpy
from numpy import zeros,matrix
from numpy import array,sinh,cosh, cos, sin, exp, arctan, linspace
from math import sqrt, pi
import pytriqs.utility.mpi as mpi
from pytriqs.gf.descriptors import Function
from triqs_ctint import SolverCore

######## parameters fixed from outside ########

U=1.0

D_continuous=False
Jperp_continuous=False

l_Ddiag=0.5
l_Doffdiag=0.5
sgn_Ddiag=1.0
sgn_Doffdiag=1.0
l_Jperp=0.5
sgn_Jperp=1.0 #screening frequency (discrete bosonic spectrum parameter, determines the sign of J_perp(tau))

w0=1.0

n_cycles = 2000

################################################
#--------------parameters fixed------------------#

mu=U/2.


beta = 20.0
n_iw = 200
delta = 0.1
n_tau_g0  = 100001
n_tau_f  = 100001
n_tau_dynamical_interactions = 1025
n_iw_dynamical_interactions = 200
n_tau_nnt = 100
n_tau_g2t = 10
n_tau_f_vertex = 61
n_tau_b_vertex = 61

# set up a few parameters

eps = 0.2

#------------------------------------------------#

N_states=1
block_names = ['up','dn']

h_int = U * n(block_names[0],0)*n(block_names[1],0)

# alpha[block][index,s]
diag = 0.5 + delta
odiag = 0.5 - delta
alpha = [ [[diag,odiag]], [[odiag,diag]] ]

gf_struct = {block_names[0] : [0], block_names[1] : [0]}

# Construct the cluster ctint solver

S = SolverCore( beta = beta,
                gf_struct = gf_struct,
                n_iw = n_iw,
                n_tau_g0  = n_tau_g0,
                n_tau_f  = n_tau_f,
                n_tau_dynamical_interactions = n_tau_dynamical_interactions,
                n_iw_dynamical_interactions = n_iw_dynamical_interactions,
                n_tau_nnt = n_tau_nnt,
                n_tau_f_vertex = n_tau_f_vertex,
                n_tau_b_vertex = n_tau_b_vertex
              )

# Initialize the Green's function

semicirc = S.G0_iw[block_names[0]].copy()
semicirc << SemiCircular(1.0)
for n,g in S.G0_iw: g << inverse(iOmega_n + mu - 1.0 * semicirc)


#------------------ set up time-dependent interactions -------------------------#
#continuous bosonic spectrum parameters
wc=1.0
s=0.2

nu = lambda n,beta : (2*n*pi/beta) #bosonic frequences
Eps=linspace(0,wc,5000) #eps grid

Dz_expr = lambda w, g : (-1.)*(.5*g**2) *((s+1)*wc**(-s-1))*sum([eps**(s+1)/(w.imag*w.imag+eps**2) for eps in Eps])*(Eps[1]-Eps[0]) #Dz away from 0
Dz0 = lambda g : (-1.)*(.5*g**2)*(s+1.)/s/wc #Dz at zero (Dz_expr ill-defined at w=0)
Dz_expr_reg = lambda w, g: Dz_expr(w,g) if abs(w)>0.0 else Dz0(g) #Dz for any w

for i in range(N_states):
  for j in range(N_states):
    if Jperp_continuous:
      S.Jperp_iw[i,j] << Function(lambda w : sgn_Jperp*4.0*Dz_expr_reg(w, l_Jperp))
    else:
      S.Jperp_iw[i,j] << sgn_Jperp*l_Jperp**2*(inverse(iOmega_n-w0)-inverse(iOmega_n+w0))

if D_continuous:
  for b1 in block_names:
    for b2 in block_names:
      for i in range(N_states):
        for j in range(N_states):
          if i==j and b1==b2:
            S.D0_iw[b1+"|"+b2][i,j] << Function(lambda w : sgn_Ddiag * Dz_expr_reg(w, l_Ddiag))
          else:
            S.D0_iw[b1+"|"+b2][i,j] << Function(lambda w : sgn_Doffdiag * Dz_expr_reg(w, l_Doffdiag))
else:
  for b1 in block_names:
    for b2 in block_names:
      for i in range(N_states):
        for j in range(N_states):
          if i==j and b1==b2:
            S.D0_iw[b1+"|"+b2][i,j]  << sgn_Ddiag*l_Ddiag**2*(inverse(iOmega_n-w0)-inverse(iOmega_n+w0))
          else:
            S.D0_iw[b1+"|"+b2][i,j]  << sgn_Doffdiag*l_Doffdiag**2*(inverse(iOmega_n-w0)-inverse(iOmega_n+w0))

if mpi.is_master_node():
  A = HDFArchive("input_quantities.h5",'w')
  A['D_iw'] = S.D0_iw
  A['G0_iw'] = S.G0_iw
  A['Jperp_tau'] = S.Jperp_tau
  A['Jperp_iw'] = S.Jperp_iw

S.solve(h_int=h_int,
        alpha = alpha,
        n_cycles = n_cycles,
        length_cycle = 50,
        n_warmup_cycles = 500,
        random_seed = 34788,
        only_sign = False,
        measure_gw = False,
        measure_ft = False,
        measure_Mt = False,
        measure_hist = False,
        measure_nnt = False,
        measure_chipmt = False,
        measure_g2t = False,
        g2t_indep = [])

gup = GfImFreq(indices = list(range(N_states)), beta = beta, n_points = n_iw, name = "upBlock")
gdn = GfImFreq(indices = list(range(N_states)), beta = beta, n_points = n_iw, name = "dnBlock")

M_iw = BlockGf(name_list = (block_names[0],block_names[1]), block_list = (gup,gdn), make_copies = True)
GfromM_iw = BlockGf(name_list = (block_names[0],block_names[1]), block_list = (gup,gdn), make_copies = True)
F_iw = BlockGf(name_list = (block_names[0],block_names[1]), block_list = (gup,gdn), make_copies = True)
GfromF_iw = BlockGf(name_list = (block_names[0],block_names[1]), block_list = (gup,gdn), make_copies = True)

for bn in block_names:
  M_iw[bn] << Fourier(S.M_tau[bn])
  F_iw[bn] << Fourier(S.F_tau[bn])
  GfromM_iw[bn] << S.G0_shift_iw[bn] + S.G0_shift_iw[bn] * M_iw[bn] * S.G0_shift_iw[bn]
  GfromF_iw[bn] <<  S.G0_shift_iw[bn] + S.G0_shift_iw[bn] * F_iw[bn]

#print out results
if mpi.is_master_node():
  A = HDFArchive("anderson.h5",'w')
  A['average_sign'] = S.average_sign
  A['G_iw'] = S.G_iw
  A['nn_tau'] = S.nn_tau
  A['chipm_tau'] = S.chipm_tau
  A['G0_iw'] = S.G0_iw
  A['G0_shift_iw'] = S.G0_shift_iw
  A['hist'] = S.histogram
  A['M_tau'] = S.M_tau
  A['M_iw'] = M_iw
  A['F_tau'] = S.F_tau
  A['F_iw'] = F_iw
  A['GfromF_iw'] = GfromF_iw
  A['GfromM_iw'] = GfromM_iw


