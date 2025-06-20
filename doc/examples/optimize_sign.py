from triqs.Base.GF_Local import *
from triqs.Base.Archive import *
from triqs.Base.Utility import MPI
from triqs.Base.Utility.myUtils import Sum
from triqs.Base.Lattice.TightBinding import *
from triqs.Base.Lattice.TightBinding import TBLattice
from triqs.Base.SumK.SumK_Discrete_From_Lattice import *
import numpy, numpy.linalg

# Parameters
x = 0.800
y = 0.126
U = 6.20
t = 1.00
tp = -0.20
Beta = 2
mu = 1.80
N_loops = 1
N_k = 8
N_Cycles = 10000

##################################################################################
# Lattice
##################################################################################

# Hoppings
hop = { ( 0, 0) : [ [  0, -t,  0, -t,-tp,  0,  0,  0,  0],
                    [ -t,  0, -t,-tp, -t,-tp,  0,  0,  0],
                    [  0, -t,  0,  0,-tp, -t,  0,  0,  0],
                    [ -t,-tp,  0,  0, -t,  0, -t,-tp,  0],
                    [-tp, -t,-tp, -t,  0, -t,-tp, -t,-tp],
                    [  0,-tp, -t,  0, -t,  0,  0,-tp, -t],
                    [  0,  0,  0, -t,-tp,  0,  0, -t,  0],
                    [  0,  0,  0,-tp, -t,-tp, -t,  0, -t],
                    [  0,  0,  0,  0,-tp, -t,  0, -t,  0] ],
        ( 1, 0) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ],
#       ( 1, 0) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
#                   [  0,  0,  0,  0,  0,  0,  0,  0,  0],
#                   [ -t,  0,  0,-tp,  0,  0,  0,  0,  0],
#                   [  0,  0,  0,  0,  0,  0,  0,  0,  0],
#                   [  0,  0,  0,  0,  0,  0,  0,  0,  0],
#                   [-tp,  0,  0, -t,  0,  0,-tp,  0,  0],
#                   [  0,  0,  0,  0,  0,  0,  0,  0,  0],
#                   [  0,  0,  0,  0,  0,  0,  0,  0,  0],
#                   [  0,  0,  0,-tp,  0,  0, -t,  0,  0] ],
        (-1, 0) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ],
        ( 0, 1) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ],
        ( 0,-1) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ],
        ( 1, 1) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ],
        ( 1,-1) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ],
        (-1, 1) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ],
        (-1,-1) : [ [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    [  0,  0,  0,  0,  0,  0,  0,  0,  0] ] }

# Lattice
L = TBLattice(Units = [(3,0,0), (0,3,0)], Hopping = hop,
              Orbital_Positions = [(0,0,0), (1,0,0), (2,0,0),
                                   (0,1,0), (1,1,0), (2,1,0),
                                   (0,2,0), (1,2,0), (2,2,0)],
              Orbital_Names = ["1","2","3","4","5","6","7","8","9"])

# Create a sum over k
SK = SumK_Discrete_From_Lattice(L, Number_Points_in_BZ = N_k)


##################################################################################
# Construct CT-INT solver
##################################################################################

from triqs.Solvers.InteractionExpansion import Solver

S = Solver(Beta = Beta,
           GFstruct = [('up',[1,2,3,4,5,6,7,8,9]), ('down',[1,2,3,4,5,6,7,8,9])],
           N_Cycles = N_Cycles,
           N_Warmup_Cycles = 5000,
           Only_Sign = True,
           Length_Cycle = 50)

# U[a,b][i,j]
S.Umatrix = [ [numpy.zeros([9,9]), U*numpy.identity(9)],
              [U*numpy.identity(9), numpy.zeros([9,9])] ]


# Prepare G0
S.Sigma << mu
S.G << SK(S.Sigma, mu = mu)
S.G0 << inverse( inverse(S.G) + S.Sigma )
G0 = S.G0.copy()


# Return average sign
def get_sign(x, f_):
  S.alpha = [ numpy.array([[x[0]],[x[1]],[x[0]],[x[1]],[x[2]],[x[1]],[x[0]],[x[1]],[x[0]]]),
              numpy.array([[x[3]],[x[4]],[x[3]],[x[4]],[x[5]],[x[4]],[x[3]],[x[4]],[x[3]]]) ]
  S.G0 << G0
  S.Solve()
  print("x = ", x[0], "  y = ", x[1], "  Av = ", S.average_sign)
  f_.write("%f %f %f\n"%(x[0],x[1],S.average_sign))
  f_.flush()
  return -abs(S.average_sign)


# Find minimum

import scipy.optimize as opt

f = open("sign.dat", 'w')
opt.fmin(get_sign, [x,x,x,y,y,y], [f])
f.close()
