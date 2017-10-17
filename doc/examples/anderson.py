from pytriqs.gf import *
from pytriqs.archive import *
from numpy import array, zeros, identity

# set up a few parameters
U = 1.0
beta = 20.0
mu = 1.3
eps = 0.2
delta = 0.35

# Interaction matrix U[block,block][index,index]
IdU = U*identity(1)
Id0 = zeros([1,1])
Umat = [[Id0,IdU],
        [IdU,Id0]]

# alpha[block][index,s]
diag = 0.5 + delta
odiag = 0.5 - delta
alpha = [[[diag,odiag]], [[odiag,diag]]]

# Construct the interaction solver
from triqs_ctint import SolverCore
S = Solver(beta = beta, block_names = ["up","down"], block_sizes = [1,1], n_iw = 200)

# Initialize the Green's function
S.G0_iw << inverse(iOmega_n + mu - inverse(iOmega_n - eps));

# Solve!
S.solve(U = Umat,
        alpha = alpha,
        beta = beta,
        n_cycles = 10000,
        length_cycle = 50,
        n_warmup_cycles = 5000,
        random_seed = 34788,
        only_sign = False,
        measure_gw = True,
        measure_ft = False,
        measure_hist = True)

# Save in archive
A = HDFArchive("anderson_py.output.h5",'w')
A["G"] = S.G_iw
