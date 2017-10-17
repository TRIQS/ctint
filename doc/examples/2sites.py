from triqs_ctint import SolverCore
from pytriqs.gf import *
from pytriqs.archive import *
import numpy as np

# Solver
S = Solver(beta = 20.0,
           block_names = ['up','down'],
           block_sizes = [2,2],
           n_iw = 300)

# Set parameters
U = 1.0
t = 0.2

# Set U and alpha
Umat = [[np.array([[0,0],[0,0]]), np.array([[U,0],[0,U]])],
        [np.array([[U,0],[0,U]]), np.array([[0,0],[0,0]])]]

alpha = [ np.array([[0.80,0.20],[0.80,0.20]]) , np.array([[0.20,0.80],[0.20,0.80]]) ]

# Prepare G0
S.G0_iw['up'][0,0] << iOmega_n + 1.3 - inverse(iOmega_n - 0.2)
S.G0_iw['up'][0,1] << t - 0.1*inverse(iOmega_n)
S.G0_iw['up'][1,0] << t - 0.1*inverse(iOmega_n)
S.G0_iw['up'][1,1] << iOmega_n + 1.3 - inverse(iOmega_n - 0.2)
S.G0_iw['down'] << S.G0_iw['up']
S.G0_iw.invert()

# Solve!
S.solve( beta = 20.0, U = Umat, alpha = alpha, n_cycles = 100000 )

# Save stuff
A = HDFArchive("twosites.output.h5",'w')
A["G_int"] = S.G_iw

