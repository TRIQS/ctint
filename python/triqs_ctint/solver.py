from solver_core  import SolverCore
from pytriqs.gf import *
from pytriqs.utility import mpi

import numpy as np
from scipy.optimize import root

np.set_printoptions(precision=4)

# === Some utility functions

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        print(arg)

# Flatten a block vector of matrices
def flatten(Sig_HF):
    return np.array([Sig_bl.flatten() for bl, Sig_bl in Sig_HF]).flatten()

# Unflatten a block vector of matrices
def unflatten(Sig_HF_flat, gf_struct):
    offset = 0
    Sig_HF = []
    for bl, indices in gf_struct:
        N = len(indices)
        Sig_HF.append([bl, Sig_HF_flat[list(range(offset, offset + N**2))].reshape((N,N))])
        offset = offset + N**2
    return Sig_HF

# === The SolverCore Wrapper

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=500, n_tau=5001, use_D=False, use_Jperp=False, n_tau_dynamical_interactions=5001, n_iw_dynamical_interactions=500):
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
        n_tau_dynamical_interactions : int, optional
               Number of tau pts for D0_tau and jperp_tau (Default 10001)
        n_iw_dynamical_interactions : int, optional
               Number of matsubara freqs for D0_iw and jperp_iw (Default 200)
        """
        if isinstance(gf_struct,dict):
            print("WARNING: gf_struct should be a list of pairs [ [str,[int,...]], ...], not a dict")
            gf_struct = list(gf_struct.items())

        # Initialise the core solver
        SolverCore.__init__(self, beta=beta, gf_struct=gf_struct, 
                            n_iw=n_iw, n_tau=n_tau, use_D=use_D, use_Jperp=use_Jperp,
                            n_tau_dynamical_interactions=n_tau_dynamical_interactions,
                            n_iw_dynamical_interactions=n_iw_dynamical_interactions)

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
                     Note that in this Python Wrapper the alpha-tensor is optional.
                     If not given, it will be constructed from the density matrix of
                     the SC Hartree Fock solution.
        """

        h_int = params_kw['h_int']
        gf_struct = self.gf_struct
         
        if 'alpha' not in params_kw:
            # --------- Determine the alpha tensor from SC Hartree Fock ----------
            mpi_print("Determine alpha-tensor from SC Hartree Fock solution")

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

                    assert(bl1 == bl2 and bl3 == bl4)

                    idx_u1 = dict(gf_struct)[bl1].index(u1)
                    idx_u2 = dict(gf_struct)[bl2].index(u2)
                    idx_u3 = dict(gf_struct)[bl3].index(u3)
                    idx_u4 = dict(gf_struct)[bl4].index(u4)
            
                    # Full Hatree Fock Solution
                    Sig_HF[bl1][idx_u2, idx_u1] += coef * G_dens[bl3][u4, u3]
                    Sig_HF[bl3][idx_u4, idx_u3] += coef * G_dens[bl1][u2, u1]

                    # # Consider cross terms for equal blocks
                    # if bl1 == bl3:
                        # Sig_HF[bl1][idx_u4, idx_u1] -= coef * G_dens[bl3][u2, u3]
                        # Sig_HF[bl3][idx_u2, idx_u3] -= coef * G_dens[bl1][u4, u1]
            
                return Sig_HF_flat - flatten(list(Sig_HF.items()))
            
            # Invoke the root finder
            Sig_HF_init = [[bl, np.zeros((len(idx_lst), len(idx_lst)))] for bl, idx_lst in gf_struct]
            root_finder = root(f, flatten(Sig_HF_init))
            
            # Now calculate alpha from the Hartree Fock solution
            delta = params_kw.pop('delta', 0.1)
            alpha = []
            densities_HF = []
            if root_finder['success']:
                Sig_HF = unflatten(root_finder['x'], gf_struct)
                mpi_print("  -- Found Sigma_HF : ")
                for bl, Sig_HF_bl in Sig_HF:
                    mpi_print("        %s : \t"%bl + str(Sig_HF_bl).replace('\n','\n           \t'))
                G_iw = self.G0_iw.copy()
                for bl, G0_bl in self.G0_iw:
                    G_iw[bl] << inverse( inverse(G0_bl) - dict(Sig_HF)[bl] )
                    dens_HF = np.diag(G_iw[bl].density()).real
                    densities_HF.append((bl, dens_HF, sum(dens_HF)))
                    if 'up' in bl:
                        alpha.append( np.array([[n_o + delta] for n_o in dens_HF ]) )
                    elif 'dn' in bl or 'do' in bl:
                        alpha.append( np.array([[n_o - delta] for n_o in dens_HF ]) )
                    else:
                        alpha.append( np.array([[n_o] for n_o in dens_HF ]) )
            else:
                mpi_print("Could not determine Hartree Fock solution, falling back to manual alpha")
                indices = gf_struct[0][1]
                for bl, G0_bl in self.G0_iw:
                    if 'up' in bl:
                        alpha.append( [[0.5 + delta] for i in indices ] )
                    elif 'dn' in bl or 'do' in bl:
                        alpha.append( [[0.5 - delta] for i in indices ] )
                    else:
                        alpha.append( [[0.5] for i in indices ] )
            
            params_kw['alpha'] = alpha
            params_kw['n_s'] = 1
            
            mpi_print("  -- HF Densities : ")
            for bl, dens, dens_tot in densities_HF:
                mpi_print("        %s : \t%s    Total: %.4f"%(bl, dens, dens_tot))

            mpi_print("  -- Alpha Tensor : ")
            for bl, alpha_bl in zip(list(dict(gf_struct).keys()), alpha):
                mpi_print("        %s : \t"%bl + str(alpha_bl).replace('\n','\n           \t'))

        # Call the core solver's solve routine
        solve_status = SolverCore.solve(self, **params_kw)

        return solve_status
