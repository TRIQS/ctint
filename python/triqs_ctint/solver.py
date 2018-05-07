from solver_core  import SolverCore
from pytriqs.gf import *
from pytriqs.utility import mpi

import numpy as np
from scipy.optimize import root

# === Some utility functions

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        print arg

# Flatten a block vector of matrices
def flatten(Sig_HF):
    return np.array([Sig_bl.flatten() for bl, Sig_bl in Sig_HF]).flatten()

# Unflatten a block vector of matrices
def unflatten(Sig_HF_flat, gf_struct):
    offset = 0
    Sig_HF = []
    for bl, indices in gf_struct:
        N = len(indices)
        Sig_HF.append([bl, Sig_HF_flat[range(offset, offset + N**2)].reshape((N,N))])
        offset = offset + N**2
    return Sig_HF

# === The SolverCore Wrapper

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau=10001, use_D=False, use_Jperp=False):
        """
        Initialise the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : list of pairs [ [str,[int,...]], ...]
                    Structure of the Green's functions. It must be a
                    list of pairs, each containing the name of the
                    Green's function block as a string and a list of integer
                    indices.
                    For example: ``[ ['up', [0, 1, 2]], ['down', [0, 1, 2]] ]``.
        n_iw : integer, optional
               Number of Matsubara frequencies used for the Green's functions.
        n_tau : integer, optional
               Number of imaginary time points used for the Green's functions.
        use_D : bool, optional
               Use dynamic density-density interaction given via S.D0_iw[bl1, bl2][i,j]
        use_Jperp : bool, optional
               Use dynamic spin-spin interaction given via S.Jperp[i,j]
        """
        if isinstance(gf_struct,dict):
            print "WARNING: gf_struct should be a list of pairs [ [str,[int,...]], ...], not a dict"
            gf_struct = list(gf_struct.iteritems())

        # Initialise the core solver
        SolverCore.__init__(self, beta=beta, gf_struct=gf_struct, 
                            n_iw=n_iw, n_tau=n_tau, use_D=use_D, use_Jperp=use_Jperp)

        self.gf_struct = gf_struct
        self.n_iw = n_iw
        self.n_tau = n_tau
        self.use_D = use_D
        self.use_Jperp = use_Jperp

    def solve(self, **params_kw):
        """
        Solve the impurity problem.

        Parameters
        ----------
        params_kw : dict {'param':value} that is passed to the core solver.
                     The only two required parameters are
                        * `h_int`: The local interaction Hamiltonian
                        * `n_cycles`: The number of Monte-Carlo cycles
                     For the other optional parameters see documentation.
                     Note that this Python Wrapper explicitly constructs the
                     alpha-tensor from the SC Hartree Fock solution
        """

        # --------- Determine the alpha tensor from Hartree Fock ----------
        mpi_print("Determine alpha-tensor from SC Hartree Fock solution")
        
        h_int = params_kw['h_int']
        gf_struct = self.gf_struct
        
        # --- Determine the self-consistent Hartree Fock solution via root finding
        
        # Find the root of this function
        def f(Sig_HF_flat):
            Sig_HF = dict(unflatten(Sig_HF_flat, gf_struct))
            G_iw = self.G0_iw.copy()
            G_dens = {}
            for bl, G0_bl in self.G0_iw:
                G_iw[bl] << inverse( inverse(G0_bl) - Sig_HF[bl] )
                G_dens[bl] = G_iw[bl].density().real
                Sig_HF[bl][:] = 0
        
            for term, coef in h_int:
                bl1, u1 = term[0][1]
                bl2, u2 = term[3][1]
                bl3, u3 = term[1][1]
                bl4, u4 = term[2][1]
        
                # Full Hatree Fock Solution
                Sig_HF[bl1][u2, u1] += coef * G_dens[bl3][u4, u3]
                Sig_HF[bl3][u4, u3] += coef * G_dens[bl1][u2, u1]
        
            return Sig_HF_flat - flatten(list(Sig_HF.iteritems()))
        
        # Invoke the root finder
        Sig_HF_init = [[bl, np.zeros((len(idx_lst), len(idx_lst)))] for bl, idx_lst in gf_struct]
        root_finder = root(f, flatten(Sig_HF_init))
        
        # --- Determine alpha from the Hartree Fock solution
        delta = 0.2
        if root_finder['success']:
            Sig_HF = unflatten(root_finder['x'], gf_struct)
            mpi_print("  -- Found Sigma_HF : ")
            for bl, Sig_HF_bl in Sig_HF: mpi_print("    " + str(bl) + " " + str(Sig_HF_bl).replace('\n',','))
            G_iw = self.G0_iw.copy()
            alpha = []
            for bl, G0_bl in self.G0_iw:
                G_iw[bl] << inverse( inverse(G0_bl) - dict(Sig_HF)[bl] )
                s = 1 if (bl == 'up') else -1
                alpha.append( [[val.real + s * delta] for val in np.diag(G_iw[bl].density()) ] )
        else:
            mpi_print("Could not determine Hartree Fock solution, falling back to manual alpha")
            indices = gf_struct[0][1]
            alpha = [ [[0.5 + delta] for i in indices ], [[0.5 - delta] for i in indices ] ]
        
        params_kw['alpha'] = alpha
        params_kw['n_s'] = 1
        
        mpi_print("  -- Alpha Tensor : " + str(alpha))

        # Call the core solver's solve routine
        solve_status = SolverCore.solve(self, **params_kw)

        return solve_status
