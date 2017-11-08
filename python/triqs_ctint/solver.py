from .solver_core import SolverCore
from triqs.gf import *
from triqs.utility import mpi

import numpy as np
from scipy.optimize import root

np.set_printoptions(precision=4)

# === Some utility functions

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        print(arg)

# === The SolverCore Wrapper

class Solver(SolverCore):

    def __init__(self, **constr_params):
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
        # Initialise the core solver
        SolverCore.__init__(self, **constr_params)

    def solve(self, **solve_params):
        """
        Solve the impurity problem.

        Parameters
        ----------
        solve_params : dict {'param':value} that is passed to the core solver.
                     The only two required parameters are
                        * `h_int`: The local interaction Hamiltonian
                        * `n_cycles`: The number of Monte-Carlo cycles
                     For the other optional parameters see documentation.
                     Note that in this Python Wrapper the alpha-tensor is optional.
                     If not given, it will be constructed automatically.
        """

        assert 'h_int' in solve_params
        assert 'n_cycles' in solve_params

        delta = solve_params.pop('delta', 0.1)
        if 'n_s' not in solve_params: solve_params['n_s'] = 1

        h_int = solve_params['h_int']
        gf_struct = self.constr_params['gf_struct']
        n_s = solve_params['n_s']
         
        if 'alpha' not in solve_params:
            # --------- Determine the alpha tensor from SC Hartree Fock ----------
            mpi_print("Determine alpha-tensor")

            if n_s != 1:
                raise Exception("Automatic alpha determination not implemented for n_s != 1")

            def sign(a):
                if a >= 0: return 1
                if a < 0: return -1

            # Prepare the inverse of G0_iw and the known high-frequency moments
            km = {}
            for bl, g_bl in self.G0_iw:
                self.G0_iw_inv[bl] << inverse(g_bl)
                km[bl] = make_zero_tail(g_bl, 2)
                km[bl][1] = np.eye(g_bl.target_shape[0])

            # The numer of terms in h_int determines the leading dimension of alpha
            n_terms = len(list(h_int))

            # Make sure that n_s parameter is compatible with automatic alpha
            n_s = params_kw.get('n_s', 1)
            if(all(any(s in bl for s in ['up', 'dn', 'do']) for bl, idxlst in gf_struct)):
                assert n_s in [1, 2], "For spinfull systems the automatic alpha determination requires n_s to be either 1 or 2"
            else:
                assert n_s == 1, "For spinless systems the automatic alpha determination requires n_s to be 1"

            # Find the root of this function
            def f(alpha_vec):
                self.eval_count = self.eval_count + 1

                # Reshape input and prepare params
                _alpha = alpha_vec.reshape(n_terms, 2, 2, n_s)
                _params = solve_params.copy()
                _params['alpha'] = _alpha
                _params.update(self.constr_params)

                # Update the shifted Green function
                self.prepare_G0_shift_iw(**_params)

                # Calculate the shifted density
                _G0_shift_dens = { bl: g_bl.density(km[bl]).real for bl, g_bl in self.G0_shift_iw }

                # Evaluate the violation of the self-consistency condition
                _res = []
                for n, (term, coeff) in enumerate(h_int):

                    bl_0, u_cdag_0 = term[0][1]
                    bl_1, u_cdag_1 = term[1][1]
                    bl_1, u_c_1 = term[2][1]
                    bl_0, u_c_0 = term[3][1]

                    idx_cdag_0 = dict(gf_struct)[bl_0].index(u_cdag_0)
                    idx_cdag_1 = dict(gf_struct)[bl_1].index(u_cdag_1)
                    idx_c_1 = dict(gf_struct)[bl_1].index(u_c_1)
                    idx_c_0 = dict(gf_struct)[bl_0].index(u_c_0)
            
                    _res.append(_G0_shift_dens[bl_0][idx_c_0,idx_cdag_0] - _alpha[n,0,0,0])
                    _res.append(_G0_shift_dens[bl_1][idx_c_1,idx_cdag_1] - _alpha[n,1,1,0])
                    # _res.append(G0_shift_dens[bl_0][idx_c_0,idx_cdag_0] + (1 if idx_c_0 == idx_cdag_0 else 0) - alpha[n,0,0,0] + sign(coeff) * delta)
                    # _res.append(G0_shift_dens[bl_1][idx_c_1,idx_cdag_1] + (1 if idx_c_1 == idx_cdag_1 else 0) - alpha[n,1,1,0] - delta)

                    if bl_0 == bl_1:
                        _res.append(_G0_shift_dens[bl_0][idx_c_1,idx_cdag_0] - _alpha[n,0,1,0])
                        _res.append(_G0_shift_dens[bl_0][idx_c_0,idx_cdag_1] - _alpha[n,1,0,0])
                        # _res.append(G0_shift_dens[bl_0][idx_c_1,idx_cdag_0] + (1 if idx_c_1 == idx_cdag_0 else 0) - alpha[n,0,1,0] - delta)
                        # _res.append(G0_shift_dens[bl_0][idx_c_0,idx_cdag_1] + (1 if idx_c_0 == idx_cdag_1 else 0) - alpha[n,1,0,0] - sign(coeff) * delta)
                    else:
                        _res.append(_alpha[n,0,1,0])
                        _res.append(_alpha[n,1,0,0])
            
                # mpi_print("Function norm " + str(np.linalg.norm(_res)))
                return np.array(_res)
            
            # Calculate initial alpha guess from G0_iw.density(km)
            alpha_init = np.zeros((n_terms, 2, 2, n_s))
            G0_dens = { bl: g_bl.density(km[bl]).real for bl, g_bl in self.G0_iw }
            for n, (term, coeff) in enumerate(h_int):

                bl_0, u_cdag_0 = term[0][1]
                bl_1, u_cdag_1 = term[1][1]
                bl_1, u_c_1 = term[2][1]
                bl_0, u_c_0 = term[3][1]

                idx_cdag_0 = dict(gf_struct)[bl_0].index(u_cdag_0)
                idx_cdag_1 = dict(gf_struct)[bl_1].index(u_cdag_1)
                idx_c_1 = dict(gf_struct)[bl_1].index(u_c_1)
                idx_c_0 = dict(gf_struct)[bl_0].index(u_c_0)
            
                alpha_init[n,0,0,0] = G0_dens[bl_0][idx_c_0,idx_cdag_0]
                alpha_init[n,1,1,0] = G0_dens[bl_1][idx_c_1,idx_cdag_1]

                if bl_0 == bl_1:
                    alpha_init[n,0,1,0] = G0_dens[bl_0][idx_c_1,idx_cdag_0]
                    alpha_init[n,1,0,0] = G0_dens[bl_0][idx_c_0,idx_cdag_1]
                else:
                    alpha_init[n,0,1,0] = 0
                    alpha_init[n,1,0,0] = 0

            # mpi_print("Init Alpha: " + str(alpha_init[...,0]))
            alpha_vec_init = alpha_init.reshape(4 * n_terms * n_s)
            alpha = alpha_init

            if mpi.is_master_node():
                # Find alpha on master_node
                self.eval_count = 0
                root_finder = root(f, alpha_vec_init, method="hybr", tol=1e-5)
                # mpi_print("Number of calls: " + str(self.eval_count))
                
                # Reshape result, Implement Fallback solution if unsuccessful
                if root_finder['success']:
                    alpha = root_finder['x'].reshape(n_terms, 2, 2, n_s)
                else:
                    mpi_print("Could not determine alpha, falling back to G0_iw.density()")

            alpha = mpi.bcast(alpha, root=0)

            #_ Introduce alpha assymetry
            for n, (term, coeff) in enumerate(h_int):
                alpha[n,0,0,0] = alpha[n,0,0,0] - sign(coeff) * delta
                alpha[n,1,1,0] = alpha[n,1,1,0] + delta
                alpha[n,0,1,0] = alpha[n,0,1,0] + delta * (abs(alpha[n,0,1,0]) > 1e-6)
                alpha[n,1,0,0] = alpha[n,1,0,0] + sign(coeff) * delta * (abs(alpha[n,1,0,0]) > 1e-6)

            mpi_print(" --- Alpha Tensor : ")
            mpi_print(str(alpha[...,0]))
            solve_params['alpha'] = alpha

        # Call the core solver's solve routine
        solve_status = SolverCore.solve(self, **solve_params)

        return solve_status
