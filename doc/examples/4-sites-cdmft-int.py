from pytriqs.gf import *
from pytriqs.lattice.tight_binding import *
from pytriqs.lattice.bz_patch import *
from pytriqs.dos.hilbert_transform import *
from pytriqs.sumk import *
from h5 import *
import pytriqs.utility.mpi as mpi
import numpy as np

# Parameters
n_loops = 1
n_k = 8
n_cycles = 50000

#x = 0.794609
#y = 0.122401
x = 0.52
y = 0.48

U = 8.00
t = -1.00
tp = -0.00
beta = 2.00
mu = 4.00

# Green's functions indices
IND = [0,1,2,3]
IND_B = ["00","01","10","11"]

# Real space
GK_iw = BlockGf(name_list = ['00-up','01-up','10-up','11-up','00-down','01-down','10-down','11-down'],
                block_list = [GfImFreq(indices = [0], beta = beta) for i in range(8)],
                make_copies = False)

##################################################################################
# Basis change (not really useful)
##################################################################################

# Basis change matrices
P = 0.5 * numpy.array([
  [ 1, 1,-1,-1],
  [ 1, 1, 1, 1],
  [ 1,-1, 1,-1],
  [ 1,-1,-1, 1]
])
Pinv = numpy.linalg.inv(P)

# Change basis
def real_to_reciprocal(green_real, green_reciprocal):
  green_reciprocal.zero()
  for i in range(4):
    for k in range(4):
      for l in range(4):
        green_reciprocal["%s-up"%(IND_B[i])][0,0] += Pinv[i,k] * green_real['up'][IND[k],IND[l]] * P[l,i]
        green_reciprocal["%s-down"%(IND_B[i])][0,0] += Pinv[i,k] * green_real['down'][IND[k],IND[l]] * P[l,i]

def reciprocal_to_real(green_reciprocal, green_real):
  green_real.zero()
  for k in range(4):
    for l in range(4):
      for i in range(4):
        green_real['up'][IND[k],IND[l]] += P[k,i] * green_reciprocal["%s-up"%(IND_B[i])][0,0] * Pinv[i,l]
        green_real['down'][IND[k],IND[l]] += P[k,i] * green_reciprocal["%s-down"%(IND_B[i])][0,0] * Pinv[i,l]

##################################################################################
# Lattice
##################################################################################

# Hoppings
hop = { ( 0, 0) : [[0,t,tp,t],[t,0,t,tp],[tp,t,0,t],[t,tp,t,0]],
        ( 1, 0) : [[0,0,0,0],[t,0,0,tp],[tp,0,0,t],[0,0,0,0]],
        (-1, 0) : [[0,t,tp,0],[0,0,0,0],[0,0,0,0],[0,tp,t,0]],
        ( 0, 1) : [[0,0,tp,t],[0,0,t,tp],[0,0,0,0],[0,0,0,0]],
        ( 0,-1) : [[0,0,0,0],[0,0,0,0],[tp,t,0,0],[t,tp,0,0]],
        ( 1, 1) : [[0,0,0,0],[0,0,0,tp],[0,0,0,0],[0,0,0,0]],
        ( 1,-1) : [[0,0,0,0],[0,0,0,0],[tp,0,0,0],[0,0,0,0]],
        (-1, 1) : [[0,0,tp,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
        (-1,-1) : [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,tp,0,0]] }

# Lattice
#BL = BravaisLattice(units = [(2,0,0), (0,2,0)], orbital_positions = [(0,1,0), (1,1,0), (1,0,0), (0,0,0)])
#TB = TightBinding(BL, hop)
TB = TBLattice(units = [(2,0,0), (0,2,0)], hopping = hop, orbital_positions = np.array([(0,1,0), (1,1,0), (1,0,0), (0,0,0)]), orbital_names = ['0','1','2','3'])

# Create a sum over k
SK = SumkDiscreteFromLattice(TB, n_points = n_k)


##################################################################################
# Construct CT-INT solver
##################################################################################

from triqs_ctint import SolverCore

S = Solver(beta = beta, gf_struct = {'up':[0,1,2,3], 'down':[0,1,2,3]})
Sigma_iw = S.G0_iw.copy()
F_iw = S.G0_iw.copy()
G_iw = S.G0_iw.copy()
G_lattice = S.G0_iw.copy()
Sigma_iw.zero()
F_iw.zero()
G_iw.zero()
G_lattice.zero()

# U[a,b][i,j]
Umat = [ [numpy.zeros([4,4]), U*numpy.identity(4)],
              [U*numpy.identity(4), numpy.zeros([4,4])] ]

# alpha[a][i,s]
alpha = [ numpy.array([[x],[x],[x],[x]]), numpy.array([[y],[y],[y],[y]]) ]


##################################################################################
# DMFT
##################################################################################

# First guess for real-space Sigma
Sigma_iw << mu

# DMFT Loop
for it in range(n_loops):

  # Construct local real-space Green's function
  G_lattice << SK(Sigma_iw, mu = mu)

  # Self-consistency
  for n,g0 in S.G0_iw:
    g0 << inverse( inverse(G_lattice[n]) + Sigma_iw[n] )

  # Run Solver
  S.solve(U = Umat,
          alpha = alpha,
          n_cycles = n_cycles,
          length_cycle = 100,
          n_warmup_cycles = 5000,
          random_seed = 34788,
          only_sign = False,
          measure_gw = False,
          measure_ft = True,
          measure_hist = True)

  # The solver provides only F(tau)
  for n,f in F_iw:
    f << Fourier(S.F_tau[n])
  G_iw << S.G0_shift_iw + S.G0_shift_iw * F_iw
  Sigma_iw << inverse(S.G0_iw) - inverse(G_iw)

  real_to_reciprocal(G_iw, GK_iw)

  # Saves
  #S.G.save("G")
  #GK.save("GK")
  if mpi.is_master_node():
    A = HDFArchive("results.h5")
    A['Sigma-%i'%it] = Sigma_iw
    A['G-%i'%it] = G_iw
    A['G0-%i'%it] = S.G0_iw
    A['GK-%i'%it] = GK_iw
    del A

