# Copyright (c) 2018--present, The Simons Foundation
# This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

from .solver_core import SolverCore
from triqs.gf import *
from triqs.utility import mpi

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
        _alpha = x.reshape(n_terms, 2, 2, 1)
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
        _res = []
        for n, (term, coeff) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self.indices_from_quartic_term(term, gf_struct)
            _res.append(_G0_shift_dens[bl0][u0p,u0] - _alpha[n,0,0,0])
            _res.append(_G0_shift_dens[bl0][u1p,u0] * (bl0 == bl1) - _alpha[n,0,1,0])
            _res.append(_G0_shift_dens[bl0][u0p,u1] * (bl0 == bl1) - _alpha[n,1,0,0])
            _res.append(_G0_shift_dens[bl1][u1p,u1] - _alpha[n,1,1,0])

        return np.array(_res)

    def convert_alpha_to_ij(self, h_int, gf_struct, n, t):
        term, coeff = list(h_int)[n]
        bl0, bl1, u0, u0p, u1, u1p = self.indices_from_quartic_term(term, gf_struct)
        if t == 0:    # alpha_00
            return [bl1, u1, u1p, -coeff]
        elif t == 1:  # alpha_01
            return [bl1, u1, u0p,  coeff] if bl0 == bl1 else None
        elif t == 2:  # alpha_10
            return [bl1, u0, u1p,  coeff] if bl0 == bl1 else None
        elif t == 3:  # alpha_11
            return [bl0, u0, u0p, -coeff]

    def jacobi(self, x, solve_params):

        # Parameters
        gf_struct = self.constr_params['gf_struct']
        h_int = solve_params['h_int']
        n_terms = len(list(h_int))

        # Update the solve_params with new alpha
        _alpha = x.reshape(n_terms, 2, 2, 1)
        _params = solve_params.copy()
        _params['alpha'] = _alpha
        _params.update(self.constr_params)

        # Prepare the inverse of G0_iw
        for bl, g_bl in self.G0_iw:
            self.G0_iw_inv[bl] << inverse(g_bl)

        # Update G0_shift_iw with new value of alpha
        self.prepare_G0_shift_iw(**_params)

        # Set trivial, diagonal derivation:
        jac = -np.eye(len(x))

        #create List with needed indices
        bl = gf_struct[0][0]
        G = self.G0_shift_iw[bl][0,0].copy()
        L  = []
        for n, (term, coeff) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self.indices_from_quartic_term(term, gf_struct)
            L.append([bl0,u0,u0p,n,0])
            if bl0 == bl1:
                L.append([bl0,u1,u0p,n,1]) # WHY, this should be u0, u1p no?
                L.append([bl1,u1,u0p,n,2])
            L.append([bl1,u1,u1p,n,3])
        for K in L:
            i = K[3]*4+K[4]
            bl = K[0]
            l = K[1] # cdag
            h = K[2] # c
            for n in range(n_terms):
                for t in range(4):
                    c = self.convert_alpha_to_ij(h_int, gf_struct, n, t)
                    if c == None:
                        continue;
                    else:
                        [blp,z1,z2,U_tilde] = c
                        j = 4*n+t
                        G.zero()
                        if bl == blp:
                            G = -self.G0_shift_iw[blp][h,z1]*self.G0_shift_iw[blp][z2,l]*U_tilde
                    jac[i,j] += np.real(G.density())
        return jac


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
                     If not given, it will be constructed from the density matrix of
                     the SC Hartree Fock solution.
        delta : float (default 0.1)
                     The value of the delta parameter used to construct the alpha tensor.
                     The larger the value of delta, the better the Monte-Carlo sign,
                     at the cost of a larger perturbation order.
                     This parameter is only used if alpha is not given explicitly.
        """

        assert 'n_cycles' in solve_params, "Solve parameter n_cycles required"

        if 'alpha' not in solve_params:

            # Parameters
            gf_struct = self.constr_params['gf_struct']
            h_int = solve_params['h_int']
            delta = solve_params.pop('delta', 0.1)
            n_s = solve_params.get('n_s', 2)
            assert n_s in [1, 2], "Solve parameter n_s has to be either 1 or 2 for automatic alpha mode"

            # --------- Determine the alpha tensor from SC Hartree Fock ----------
            mpi_print("Determine alpha-tensor")

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
                        alpha_init[n,0,0,0] -= - sign(coeff) * delta
                        alpha_init[n,1,1,0] -= delta
                        alpha_init[n,0,1,0] -= delta * (abs(alpha_init[n,0,1,0]) > 1e-6)
                        alpha_init[n,1,0,0] -= sign(coeff) * delta * (abs(alpha_init[n,1,0,0]) > 1e-6)
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
                    alpha_init[n,0,0,0] = G0_dens[bl0][u0p,u0]
                    alpha_init[n,0,1,0] = G0_dens[bl0][u1p,u0] * (bl0 == bl1)
                    alpha_init[n,1,0,0] = G0_dens[bl0][u0p,u1] * (bl0 == bl1)
                    alpha_init[n,1,1,0] = G0_dens[bl1][u1p,u1]

            # mpi_print("Init Alpha: " + str(alpha_init[...,0]))
            alpha_vec_init = alpha_init.reshape(4 * n_terms)
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
                        alpha_sc = root_finder['x'].reshape(n_terms, 2, 2, 1)
                        found = True
                    elif count > 100:
                        mpi_print("Could not determine alpha, falling back to G0_iw.density()")
                        alpha_sc = alpha_init
                        found = True
                    else:
                        from numpy import random
                        mpi_print("Could not determine alpha, Try again step %s" % count)
                        for i in range(len(alpha_vec_init)):
                            alpha_vec_init[i] = alpha_vec_init[i] + (-0.5+random.rand())
                #_ Introduce alpha assymetry
                for n, (term, coeff) in enumerate(h_int):
                    for s in range(n_s):
                        alpha[n,0,0,s] = alpha_sc[n,0,0,0] - sign(coeff) * delta * (1 - 2*s)
                        alpha[n,1,1,s] = alpha_sc[n,1,1,0] + delta * (1 - 2*s)
                        alpha[n,0,1,s] = alpha_sc[n,0,1,0] + delta * (1 - 2*s) * (abs(alpha_sc[n,0,1,0]) > 1e-6)
                        alpha[n,1,0,s] = alpha_sc[n,1,0,0] + sign(coeff) * delta * (1 - 2*s) * (abs(alpha_sc[n,1,0,0]) > 1e-6)

            alpha = mpi.bcast(alpha, root=0)

            # Make sure to set n_s as provided by the user
            solve_params['n_s'] = n_s

            mpi_print(" --- Alpha Tensor : ")
            if n_s == 1:
                mpi_print(str(alpha[...,0]))
            else:
                mpi_print("Alpha Tensor s = 1:")
                mpi_print(str(alpha[...,0]))
                mpi_print("Alpha Tensor s = 2:")
                mpi_print(str(alpha[...,1]))
            solve_params['alpha'] = alpha

        solve_status = SolverCore.solve(self, **solve_params)
        return solve_status
