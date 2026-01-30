# Copyright (c) 2018--present, The Simons Foundation
# This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

from .solver_core import SolverCore
from triqs.gf import *
from triqs.utility import mpi
from triqs_hartree_fock import ImpuritySolver as HFSolver

import numpy as np


# === Some utility functions

# print on master node
def mpi_print(arg):
    np.set_printoptions(precision=4)
    if mpi.is_master_node():
        po = np.get_printoptions()
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


    def _indices_from_quartic_term(self, term):
        """
        Extract indices from a quartic operator term.

        For a term: U_l * cdag_[bl0,u_0] cdag_[bl1,u_1] c_[bl1,u_1p] c_[bl0,u_0p]
        Returns: [bl0, bl1, u0, u0p, u1, u1p]
        """
        bl0, u0 = term[0][1]
        bl1, u1 = term[1][1]
        bl1p, u1p = term[2][1]
        bl0p, u0p = term[3][1]
        assert bl0 == bl0p and bl1 == bl1p
        return [bl0, bl1, u0, u0p, u1, u1p]

    def find_alpha_from_HF_solver(self, solve_params):
        """
        Determine the alpha tensor using triqs_hartree_fock.ImpuritySolver.

        The HF solver finds self-consistent G_iw from which we extract
        the density matrix elements needed for the alpha tensor.
        """
        mpi_print("Determine alpha-tensor using triqs_hartree_fock")

        gf_struct = self.constr_params['gf_struct']
        beta = self.constr_params['beta']
        n_iw = self.constr_params['n_iw']
        h_int = solve_params['h_int']
        delta = solve_params.pop('delta', [0.1, 0.1])
        n_s = solve_params.get('n_s', 2)
        assert n_s in [1, 2], "Solve parameter n_s has to be either 1 or 2 for automatic alpha mode"

        def sign(a):
            return 1 if a >= 0 else -1

        # The number of terms in h_int determines the leading dimension of alpha
        n_terms = len(list(h_int))

        # Create HF solver instance
        hf_solver = HFSolver(
            gf_struct=gf_struct,
            beta=beta,
            n_iw=n_iw,
            dc=False,
            force_real=True
        )

        # Copy G0_iw to the HF solver
        hf_solver.G0_iw << self.G0_iw

        # Initialize Sigma_HF from previous alpha if available
        if self.last_solve_params is not None:
            mpi_print("Initializing HF solver from previous iteration")
            alpha_prev = self.last_solve_params['alpha']
            self._initialize_hf_sigma_from_alpha(hf_solver, h_int, alpha_prev, delta)

        # Run self-consistent HF (only on master node, HF solver handles MPI internally)
        hf_solver.solve(
            h_int=h_int,
            with_fock=True,
            one_shot=False,
            method='krylov',
            tol=1e-10
        )

        # Map HF density to alpha tensor (self-consistent values without delta shift)
        alpha_sc = np.empty((n_terms, 3))
        for n, (term, coeff) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self._indices_from_quartic_term(term)
            alpha_sc[n, 0] = hf_solver.density[bl0][u0p, u0]
            alpha_sc[n, 1] = hf_solver.density[bl0][u0p, u1] * (bl0 == bl1)
            alpha_sc[n, 2] = hf_solver.density[bl1][u1p, u1]

        alpha_sc = mpi.bcast(alpha_sc, root=0)

        # Apply delta shift for n_s spin components
        alpha = np.zeros((n_terms, 2, 2, n_s))
        for n, (term, coeff) in enumerate(h_int):
            for _s in range(n_s):
                s = 1 - 2 * _s
                alpha[n, 0, 0, _s] = alpha_sc[n, 0] - sign(coeff) * s * delta[0]
                alpha[n, 0, 1, _s] = alpha_sc[n, 1] + s * delta[1] * (abs(alpha_sc[n, 1]) > 1e-6)
                alpha[n, 1, 0, _s] = alpha_sc[n, 1] + sign(coeff) * s * delta[1] * (abs(alpha_sc[n, 1]) > 1e-6)
                alpha[n, 1, 1, _s] = alpha_sc[n, 2] + s * delta[0]

        # Broadcast result (HF solver may have already done this, but ensure consistency)
        alpha = mpi.bcast(alpha, root=0)

        # Make sure to set n_s as provided by the user
        solve_params['n_s'] = n_s

        return alpha

    def _initialize_hf_sigma_from_alpha(self, hf_solver, h_int, alpha_prev, delta):
        """
        Initialize HF solver's Sigma_HF from a previous alpha tensor.

        This provides a warm start for the self-consistency loop by
        converting alpha back to an approximate self-energy.
        """
        gf_struct = self.constr_params['gf_struct']

        def sign(a):
            return 1 if a >= 0 else -1

        # Undo delta shift to get self-consistent alpha
        if alpha_prev.shape[-1] > 1:
            alpha_sc = np.mean(alpha_prev, axis=-1)
        else:
            alpha_sc = alpha_prev[..., 0].copy()
            # Undo the delta shift for n_s=1 case
            for n, (term, coeff) in enumerate(h_int):
                alpha_sc[n, 0, 0] -= -sign(coeff) * delta[0]
                alpha_sc[n, 0, 1] -= delta[1] * (abs(alpha_sc[n, 0, 1] - delta[1]) > 1e-6)
                alpha_sc[n, 1, 0] -= sign(coeff) * delta[1] * (abs(alpha_sc[n, 1, 0] - delta[1]) > 1e-6)
                alpha_sc[n, 1, 1] -= delta[0]

        # Reset Sigma_HF to zero
        for bl, _ in gf_struct:
            hf_solver.Sigma_HF[bl][:] = 0.0

        # Build approximate Sigma_HF from alpha
        # The mapping follows the Hartree-Fock equations
        for n, (term, coeff) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self._indices_from_quartic_term(term)

            # Hartree terms: Sigma[bl0][u0,u0p] += coeff * alpha[n,1,1] (density of other pair)
            #                Sigma[bl1][u1,u1p] += coeff * alpha[n,0,0]
            hf_solver.Sigma_HF[bl0][u0p, u0] += coeff * alpha_sc[n, 1, 1]
            hf_solver.Sigma_HF[bl1][u1p, u1] += coeff * alpha_sc[n, 0, 0]

            # Fock terms (if same block)
            if bl0 == bl1:
                hf_solver.Sigma_HF[bl0][u0p, u1] -= coeff * alpha_sc[n, 1, 0]
                hf_solver.Sigma_HF[bl0][u1p, u0] -= coeff * alpha_sc[n, 0, 1]

    def trivial_alpha(self, solve_params):
        h_int = solve_params['h_int']
        n_terms = len(list(h_int))
        delta = solve_params.pop('delta', [0.5 + 1e-2, 1e-2])

        assert solve_params['n_s'] == 2
        alpha = np.zeros((n_terms, 2, 2, 2))
        for l, (term, _) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self._indices_from_quartic_term(term)

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
            delta = solve_params.get('delta', [0.1, 0.1])
            try:
                iter(delta) # check if delta is iterable
                assert len(delta) == 2, "delta can only have two components"
            except TypeError:
                # catch the non-iterable case and convert to list
                solve_params['delta'] = [delta, delta]

            alpha_mode = solve_params.pop('alpha_mode', "automatic")
            if alpha_mode == "automatic":
                alpha = self.find_alpha_from_HF_solver(solve_params)
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
