from .solver_core import SolverCore
from triqs.gf import *
from triqs.utility import mpi

import itertools
import numpy as np
from scipy.optimize import root


# === Some utility functions

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        po = np.get_printoptions()
        np.set_printoptions(precision=4)
        print(arg)
        np.set_printoptions(**po)

# === The SolverCore Wrapper

class Solver(SolverCore):

    def __init__(self, **constr_params):
        """
        Initialise the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : list of pairs [ (str,int), ...]
                    Structure of the Green's functions. It must be a
                    list of pairs, each containing the name of the
                    Green's function block as a string and the size of that block.
                    For example: ``[ ('up', 3), ('down', 3) ]``.
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
        constr_params['gf_struct'] = fix_gf_struct_type(constr_params['gf_struct'])

        # Initialise the core solver
        SolverCore.__init__(self, **constr_params)


    # Helper Function
    #
    # For a given quartic operator
    #
    #   U_l * cdag_[bl0,u_0] cdag_[bl1,u_1] c_[bl1,u_1p] c_[bl0,u_0p]
    #
    # return the list of indicies
    #
    #   [bl0, bl1, u0, u0p, u1, u1p]
    #
    def indices_from_quartic_term(self, term, gf_struct):
        bl0, u0 = term[0][1]
        bl1, u1 = term[1][1]
        bl1p, u1p = term[2][1]
        bl0p, u0p = term[3][1]
        assert bl0 == bl0p and bl1 == bl1p
        return [bl0, bl1, u0, u0p, u1, u1p]

    def f(self, x, solve_params):

        # Parameters
        gf_struct = self.constr_params['gf_struct']
        h_int = solve_params['h_int']
        n_terms = len(list(h_int))

        # Update the solve_params with new alpha
        tmp = x.reshape(n_terms, 3)
        _alpha = np.empty((n_terms, 2, 2, 1))
        _alpha[:,0,0,0] = tmp[:,0]
        _alpha[:,0,1,0] = tmp[:,1]
        _alpha[:,1,0,0] = tmp[:,1]
        _alpha[:,1,1,0] = tmp[:,2]
        _params = solve_params.copy()
        _params['alpha'] = _alpha
        _params.update(self.constr_params)

        # Prepare the inverse of G0_iw and the known high-frequency moments
        km = {}
        for bl, g_bl in self.G0_iw:
            self.G0_iw_inv[bl] << inverse(g_bl)
            km[bl] = make_zero_tail(g_bl, 2)
            km[bl][1] = np.eye(g_bl.target_shape[0])

        # Update G0_shift_iw with new value of alpha
        self.prepare_G0_shift_iw(**_params)

        # Precalculate the G0_shift_iw densities
        _G0_shift_dens = { bl: g_bl.density(km[bl]).real for bl, g_bl in self.G0_shift_iw }

        # Evaluate the violation of the self-consistency condition
        _res = np.empty((n_terms,3))
        for n, (term, coeff) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self.indices_from_quartic_term(term, gf_struct)
            _res[n,0] = _G0_shift_dens[bl0][u0p,u0] - _alpha[n,0,0,0]
            #_res[n,1] = _G0_shift_dens[bl0][u1p,u0] * (bl0 == bl1) - _alpha[n,0,1,0]
            _res[n,1] = _G0_shift_dens[bl0][u0p,u1] * (bl0 == bl1) - _alpha[n,1,0,0]
            _res[n,2] = _G0_shift_dens[bl1][u1p,u1] - _alpha[n,1,1,0]

        return _res.reshape(n_terms * 3)

    def jacobi(self, x, solve_params):

        # Parameters
        gf_struct = self.constr_params['gf_struct']
        h_int = solve_params['h_int']
        n_terms = len(list(h_int))

        # Update the solve_params with new alpha
        tmp = x.reshape(n_terms, 3)
        _alpha = np.empty((n_terms, 2, 2, 1))
        _alpha[:,0,0,0] = tmp[:,0]
        _alpha[:,0,1,0] = tmp[:,1]
        _alpha[:,1,0,0] = tmp[:,1]
        _alpha[:,1,1,0] = tmp[:,2]
        _params = solve_params.copy()
        _params['alpha'] = _alpha
        _params.update(self.constr_params)

        # Prepare the inverse of G0_iw
        for bl, g_bl in self.G0_iw:
            self.G0_iw_inv[bl] << inverse(g_bl)

        # Update G0_shift_iw with new value of alpha
        self.prepare_G0_shift_iw(**_params)
        assert is_gf_hermitian(self.G0_shift_iw)

        # Initialize the Jacobi matrix
        jac = np.zeros((n_terms * 3, n_terms * 3))

        #create List with needed indices
        bl = gf_struct[0][0]
        G = self.G0_shift_iw[bl][0,0].copy()

        delta = lambda a, b: float(a == b)
        for (l, (term, coeff)), a, b in itertools.product(enumerate(h_int), range(2), range(2)):
            for (l_, (term_, coeff_)), a_, b_ in itertools.product(enumerate(h_int), range(2), range(2)):
                if (a == 0 and b == 1):
                    continue
                i = 3 * l + a + b
                j = 3 * l_ + a_ + b_

                bl0, bl1, u0, u0p, u1, u1p = self.indices_from_quartic_term(term, gf_struct)
                bl0_, bl1_, u0_, u0p_, u1_, u1p_ = self.indices_from_quartic_term(term_, gf_struct)
                bl, bl_ = [bl0, bl1], [bl1_, bl0_]
                u, u_ = [u0, u1], [u1_, u0_]
                up, up_ = [u0p, u1p], [u1p_, u0p_]

                dens = 0
                if bl[a] == bl[b] and bl[a] == bl_[b_] and bl_[a_] == bl[b]:
                    pm = 2.*(delta(a_, b_) - .5)
                    G << self.G0_shift_iw[bl[b]][up[b],u_[a_]] * self.G0_shift_iw[bl[a]][up_[b_],u[a]]
                    dens = pm * coeff_ * np.real(G.density())
                    jac[i,j] += dens
                jac[i,j] += -delta(l, l_) * delta(a, a_) * delta(b, b_)
        return jac

    def find_alpha_from_self_consistent_HF(self, solve_params):
        # --------- Determine the alpha tensor from SC Hartree Fock ----------
        mpi_print("Determine alpha-tensor")
        assert not (self.constr_params['use_D'] or self.constr_params['use_Jperp']), \
            "Automatic determination of alpha tensor does not work with D0_iw or Jperp_iw"

        gf_struct = self.constr_params['gf_struct']
        h_int = solve_params['h_int']
        delta = solve_params.pop('delta', [0.1, 0.1])
        n_s = solve_params.get('n_s', 1)
        assert n_s in [1, 2], "Solve parameter n_s has to be either 1 or 2 for automatic alpha mode"

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

        # We always solve the self-consistency assuming n_s == 1
        # If n_s > 1 we use this only in the post-processing of alpha
        solve_params['n_s'] = 1
        if not self.last_solve_params is None:
            mpi_print("Reusing alpha from previous iteration")
            alpha_init = self.last_solve_params['alpha']
            if alpha_init.shape[-1] == 1:
                # undo the delta shift
                for n, (term, coeff) in enumerate(h_int):
                    alpha_init[n,0,0,0] -= - sign(coeff) * delta[0]
                    alpha_init[n,0,1,0] -= delta[1] * (abs(alpha_init[n,0,1,0] - delta[1]) > 1e-6)
                    alpha_init[n,1,0,0] -= sign(coeff) * delta[1] * (abs(alpha_init[n,1,0,0] - delta[1]) > 1e-6)
                    alpha_init[n,1,1,0] -= delta[0]
            else:
                # Project the supplied alpha on n_s = 1 in case the provided one was n_s = 2
                # This will naturally remove the opposite delta
                alpha_init = np.mean(alpha_init, axis=-1).reshape(n_terms, 2, 2, 1)
        else:
            # Calculate initial alpha guess from G0_iw.density(km)
            alpha_init = np.zeros((n_terms, 2, 2, 1))
            # Precalculate the G0_iw densities
            G0_dens = { bl: g_bl.density(km[bl]).real for bl, g_bl in self.G0_iw }
            for n, (term, coeff) in enumerate(h_int):
                bl0, bl1, u0, u0p, u1, u1p = self.indices_from_quartic_term(term, gf_struct)
                alpha_init[n,0,0,0] = G0_dens[bl0][u1p,u1]
                alpha_init[n,0,1,0] = G0_dens[bl0][u0p,u1] * (bl0 == bl1)
                alpha_init[n,1,0,0] = G0_dens[bl0][u1p,u0] * (bl0 == bl1)
                alpha_init[n,1,1,0] = G0_dens[bl1][u0p,u0]

        # mpi_print("Init Alpha: " + str(alpha_init[...,0]))
        alpha_vec_init = np.empty((n_terms, 3))
        alpha_vec_init[:,0] = alpha_init[:,0,0,0]
        alpha_vec_init[:,1] = alpha_init[:,1,0,0]
        alpha_vec_init[:,2] = alpha_init[:,1,1,0]
        alpha_vec_init = alpha_vec_init.reshape(n_terms * 3)
        alpha = np.zeros((n_terms, 2, 2, n_s))

        if mpi.is_master_node():
            found = False
            count = 0
            while (found == False):
                count  +=1

                # Find alpha on master_node
                if solve_params.pop('use_jacobi', True):
                    mpi_print("Using Jacobi-Matrix for root search")
                    root_finder = root(self.f, alpha_vec_init,args=(solve_params),jac=self.jacobi,method="hybr")
                else:
                    root_finder = root(self.f, alpha_vec_init,args=(solve_params),method="hybr")
                # Reshape result, Implement Fallback solution if unsuccessful
                if root_finder['success']:
                    alpha_sc = root_finder['x'].reshape(n_terms, 3)
                    found = True
                elif count > 100:
                    mpi_print("Could not determine alpha, falling back to G0_iw.density()")
                    alpha_sc = alpha_vec_init
                    found = True
                else:
                    from numpy import random
                    mpi_print("Could not determine alpha, Try again step %s" % count)
                    for i in range(len(alpha_vec_init)):
                        alpha_vec_init[i] = alpha_vec_init[i] + (-0.5+random.rand())
            #_ Introduce alpha assymetry
            for n, (term, coeff) in enumerate(h_int):
                for _s in range(n_s):
                    s = 1 - 2 * _s
                    alpha[n,0,0,_s] = alpha_sc[n,0] - sign(coeff) * s * delta[0]
                    alpha[n,0,1,_s] = alpha_sc[n,1] + s * delta[1] * (abs(alpha_sc[n,1]) > 1e-6)
                    alpha[n,1,0,_s] = alpha_sc[n,1] + sign(coeff) * s * delta[1] * (abs(alpha_sc[n,1]) > 1e-6)
                    alpha[n,1,1,_s] = alpha_sc[n,2] + s * delta[0]

        alpha = mpi.bcast(alpha, root=0)

        # Make sure to set n_s as provided by the user
        solve_params['n_s'] = n_s

        return alpha

    def trivial_alpha(self, solve_params):
        gf_struct = self.constr_params['gf_struct']
        use_D = self.constr_params['use_D']
        use_Jperp = self.constr_params['use_Jperp']
        h_int = solve_params['h_int']
        n_terms = len(list(h_int))
        delta = solve_params.pop('delta', [0.5 + 1e-2, 1e-2])

        assert solve_params['n_s'] == 2
        alpha = np.zeros((n_terms + int(use_D) + int(use_Jperp), 2, 2, 2))
        for l, (term, coeff) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self.indices_from_quartic_term(term, gf_struct)

            # on-site density-density
            if bl0 != bl1 and u0 == u1 and u0p == u1p and u0 == u0p and u1 == u1p:
                alpha_s = lambda s: np.array([[ 0.5 + s*delta[0], 0.0              ],
                                              [ 0.0             , 0.5 - s*delta[0] ]])
                alpha[l,...,0] = alpha_s(+1)
                alpha[l,...,1] = alpha_s(-1)
            # inter-site density-density "up-down"
            elif bl0 != bl1 and u0 != u1 and u0p != u1p and u0 == u0p and u1 == u1p:
                alpha_s = lambda s: np.array([[ 0.5 + s*delta[0], 0.0              ],
                                              [ 0.0             , 0.5 - s*delta[0] ]])
                alpha[l,...,0] = alpha_s(+1)
                alpha[l,...,1] = alpha_s(-1)
            # inter-site density-density "up-up" or "down-down"
            elif bl0 == bl1 and u0 != u1 and u0p != u1p and u0 == u0p and u1 == u1p:
                alpha_s = lambda s: np.array([[ 0.5 + s*delta[0],       s*delta[1] ],
                                              [     - s*delta[1], 0.5 - s*delta[0] ]])
                alpha[l,...,0] = alpha_s(+1)
                alpha[l,...,1] = alpha_s(-1)
            # spin-flip
            elif bl0 != bl1 and u0 != u1 and u0p != u1p and u0 == u1p and u1 == u0p:
                assert False, "Spin-flip terms are not yet treated"
            # pair-hopping
            elif bl0 != bl1 and u0 == u1 and u0p == u1p and u0 != u0p and u1 != u1p:
                assert False, "Pair-hopping terms are not yet treated"
            else:
                assert False, "I don't know this type of term"

        if use_D:
            # F.F. Assaad, T.C. Lang, Phys. Rev. B 76, 035116 (2007) -- Eq. (36)
            alpha_s = lambda s: np.array([[ 0.5 + s*delta[0], 0.0              ],
                                          [ 0.0             , 0.5 + s*delta[0] ]])
            alpha[n_terms,...,0] = alpha_s(+1)
            alpha[n_terms,...,1] = alpha_s(-1)

        if use_Jperp:
            alpha[n_terms + int(use_D),...,0] = 0.0
            alpha[n_terms + int(use_D),...,1] = 0.0

        return alpha

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

        assert 'n_cycles' in solve_params, "Solve parameter n_cycles required"

        if 'alpha' not in solve_params:

            # Parameters
            delta = solve_params.get('delta', [0.1, 0.1])
            try:
                iter(delta) # check if delta is iterable
                assert len(delta) == 2, "delta can only have two components"
            except TypeError:
                # catch the non-iterable case and convert to list
                solve_params['delta'] = [delta, delta]

            alpha_mode = solve_params.pop('alpha_mode', "automatic")
            if alpha_mode == "automatic":
                alpha = self.find_alpha_from_self_consistent_HF(solve_params)
            elif alpha_mode == "trivial":
                alpha = self.trivial_alpha(solve_params)
            else:
                assert False, f"No such alpha_mode: {alpha_mode}"

            mpi_print(" --- Alpha Tensor : ")
            if solve_params['n_s'] == 1:
                mpi_print(str(alpha[...,0]))
            else:
                mpi_print("Alpha Tensor s = 1:")
                mpi_print(str(alpha[...,0]))
                mpi_print("Alpha Tensor s = 2:")
                mpi_print(str(alpha[...,1]))
            solve_params['alpha'] = alpha

        solve_status = SolverCore.solve(self, **solve_params)
        return solve_status
