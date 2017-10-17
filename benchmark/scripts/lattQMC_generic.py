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

################################################
# Hamiltonian creator
def initCubicTBH(Nx, Ny, Nz, eps, t, cyclic=True):
  H = [[0 for j in range(Nx*Ny*Nz)] for i in range(Nx*Ny*Nz)]  
  for i in range(Nx*Ny*Nz):
    H[i][i]=eps    
  for i in range(Nx):
    for j in range(Ny):
      for k in range(Nz): 
        if Nx>1:
          if i+1==Nx:
            if cyclic: H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + j*Nx ] = t
          else:
            H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + j*Nx + i+1 ] = t
        
          if i==0:
            if cyclic: H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + j*Nx + Nx-1] = t
          else:
            H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + j*Nx + i-1] = t  
            
        if Ny>1:
          if j+1==Ny:
            if cyclic: H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + i ] = t
          else:  
            H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + (j+1)*Nx + i ] = t
        
          if j==0:
            if cyclic: H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + (Ny-1)*Nx + i ] = t
          else:
            H [ k*Nx*Ny+j*Nx+i ]  [ k*Nx*Ny + (j-1)*Nx + i ] = t

        if Nz>1:
          if (k+1==Nz): 
            if cyclic: H [ k*Nx*Ny+j*Nx+i ]  [ j*Nx + i ] = t
          else:
            H [ k*Nx*Ny+j*Nx+i ]  [ (k+1)*Nx*Ny + j*Nx + i ] = t 
            
          if k==0:
            if cyclic: H [ k*Nx*Ny+j*Nx+i ]  [ (Nz-1)*Nx*Ny + j*Nx + i ] = t
          else:
            H [ k*Nx*Ny+j*Nx+i ]  [ (k-1)*Nx*Ny + j*Nx + i ] = t
    
  return H 

def delta(i,j):
  if i==j: return 1
  return 0

################################################


######## parameters fixed from outside ########

U=a
mu=a
beta=a

n_cycles = 10000

################################################

#------------------------------------------------#
# PREPARE H_0, H_int 

Nx=Ny=a
N_states=Nx*Ny

block_names = ['up','dn']
gf_struct = {block_names[0] : range(N_states), block_names[1] : range(N_states)}

h_int = U * n(block_names[0],0)*n(block_names[1],0)
for i in range(1,N_states):
  h_int += U * n(block_names[0],i)*n(block_names[1],i)

eps = -mu
t = -0.25
cyclic=a
H = initCubicTBH(Nx, Ny, 1, eps, t, cyclic)

#------------------------------------------------#

# Construct the cluster ctint solver
n_iw = 200
n_tau = 10001

S = SolverCore( beta = beta, 
                gf_struct = gf_struct,
                n_iw = n_iw,  
                n_tau_g0  = n_tau,
                n_tau_f  = n_tau,
                n_tau_dynamical_interactions = 3,
                n_iw_dynamical_interactions = 1,
                n_tau_nnt = 500,
                n_tau_f_vertex = 10,
                n_tau_b_vertex = 10
              )

# Initialize the Green's function
for b in block_names:
  for i in range(N_states):
    for j in range(N_states):    
      S.G0_iw[b][i,j] = delta(i,j)*iOmega_n - H[i][j]
  S.G0_iw[b].invert()

if mpi.is_master_node():
  A = HDFArchive("input_quantities.h5",'w')
  A['G0_iw'] = S.G0_iw

# alpha[block][index,s]
N_s = 2
delta = 0.1
alpha = [ [ [ 0.5 + delta*(-1)**(s+sig) for s in range(N_s)] for i in range(N_states)] for sig in range(2) ]

g2t_indep = []

# Solve!
S.solve(h_int=h_int,
        alpha = alpha,
        n_cycles = n_cycles,
        length_cycle = 50,
        n_warmup_cycles = 5000,
        random_seed = 34788,
        only_sign = False,
        measure_sign = True, 
        measure_gw = False,
        measure_ft = False,
        measure_Mt = False,
        measure_hist = True,
        measure_nnt = False,
        measure_nn = True,
        measure_chipmt = False,
        measure_g2t = False,
        g2t_indep = g2t_indep)

gup = GfImFreq(indices = range(N_states), beta = beta, n_points = n_iw, name = "upBlock")
gdn = GfImFreq(indices = range(N_states), beta = beta, n_points = n_iw, name = "dnBlock")

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
  A['hist'] = S.histogram
  A['nn'] = S.nn
#  A['G_iw'] = S.G_iw
#  A['nn_tau'] = S.nn_tau
#  A['chipm_tau'] = S.chipm_tau
#  A['G0_iw'] = S.G0_iw
#  A['G0_shift_iw'] = S.G0_shift_iw
#  A['M_tau'] = S.M_tau
#  A['M_iw'] = M_iw
#  A['F_tau'] = S.F_tau
#  A['F_iw'] = F_iw
#  A['GfromF_iw'] = GfromF_iw
#  A['GfromM_iw'] = GfromM_iw
