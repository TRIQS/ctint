# Copyright (c) 2018--present, The Simons Foundation
# This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

from .solver_core import SolverCore
from .version import gtau_is_complex, interaction_is_complex
from triqs.gf import *
from triqs.utility import mpi
from triqs_hartree_fock import ImpuritySolver as HFSolver

import numpy as np


# === Some utility functions

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        with np.printoptions(precision=4):
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
        gf_struct : list of pairs [ (str,int), ...]
                    Structure of the Green's functions. It must be a
                    list of pairs, each containing the name of the
                    Green's function block as a string and the size of that block.
                    For example: ``[ ('up', 3), ('down', 3) ]``.
        dlr_wmax : float
               DLR bandwidth cutoff w_max (= Lambda / beta).
        dlr_eps : float, optional
               DLR error tolerance epsilon. Default: 1e-10
        n_tau : integer, optional
               Number of imaginary time points used for the Green's functions.
        use_D : bool, optional
               Use dynamic density-density interaction given via S.D0_iw[bl1, bl2][i,j]
        use_Jperp : bool, optional
               Use dynamic spin-spin interaction given via S.Jperp[i,j]
        n_tau_dynamical_interactions : int, optional
               Number of tau pts for D0_tau and jperp_tau (Default 10001)
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
        return (bl0, bl1, u0, u0p, u1, u1p)

    def find_alpha_from_HF_solver(self, solve_params):
        """
        Determine the alpha tensor using triqs_hartree_fock.ImpuritySolver.

        The HF solver finds self-consistent G_iw from which we extract
        the density matrix elements needed for the alpha tensor.
        """
        mpi_print("Determine alpha-tensor using triqs_hartree_fock")

        gf_struct = self.constr_params['gf_struct']
        beta = self.constr_params['beta']
        h_int = solve_params['h_int']
        delta = solve_params.pop('delta', [0.1, 0.1])
        n_s = solve_params.get('n_s', 2)
        assert n_s in [1, 2], "Solve parameter n_s has to be either 1 or 2 for automatic alpha mode"

        # The number of terms in h_int determines the leading dimension of alpha
        n_terms = len(list(h_int))

        # Create HF solver instance with same DLR parameters as ctint
        dlr_wmax = self.constr_params['dlr_wmax']
        dlr_eps = self.constr_params['dlr_eps']
        hf_solver = HFSolver(
            gf_struct=gf_struct,
            beta=beta,
            w_max=dlr_wmax,
            eps=dlr_eps,
            dc=False,
            force_real=not gtau_is_complex
        )

        # Copy G0_iw to the HF solver (same DLR mesh)
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

        # Determine number of D0 alpha entries
        n_D0_total = 0
        if self.constr_params['use_D']:
            block_names = [bl for bl, _ in gf_struct]
            n_bl = len(block_names)
            R = gf_struct[0][1]
            n_D0_total = n_bl * n_bl * R * R

        # Build alpha tensor from HF density with delta shift for n_s spin components
        alpha = np.zeros((n_terms + n_D0_total, 2, 2, n_s), dtype=complex if gtau_is_complex else float)
        for n, (term, coeff) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self._indices_from_quartic_term(term)
            n00 = hf_solver.density[bl0][u0p, u0]
            n01 = hf_solver.density[bl0][u0p, u1] * (bl0 == bl1)
            n11 = hf_solver.density[bl1][u1p, u1]
            has_offdiag = abs(n01) > 1e-6
            for s in range(n_s):
                sgn = 1 - 2 * s # delta sign for each aux spin component
                alpha[n, 0, 0, s] = n00 - np.sign(coeff) * sgn * delta[0]
                alpha[n, 0, 1, s] = n01 + sgn * delta[1] * has_offdiag
                alpha[n, 1, 0, s] = n01 + np.sign(coeff) * sgn * delta[1] * has_offdiag
                alpha[n, 1, 1, s] = n11 + sgn * delta[0]

        # Fill D0 alpha entries with per-orbital diagonal densities
        if n_D0_total > 0:
            for ibl1, (bl1, _) in enumerate(gf_struct):
                for ibl2, (bl2, _) in enumerate(gf_struct):
                    for a in range(R):
                        for b in range(R):
                            d = ibl1 * n_bl * R * R + ibl2 * R * R + a * R + b
                            for s in range(n_s):
                                sgn = 1 - 2 * s
                                alpha[n_terms + d, 0, 0, s] = hf_solver.density[bl1][a, a] + sgn * delta[0]
                                alpha[n_terms + d, 1, 1, s] = hf_solver.density[bl2][b, b] + sgn * delta[0]

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

        # Undo delta shift to get self-consistent alpha
        if alpha_prev.shape[-1] > 1:
            alpha_sc = np.mean(alpha_prev, axis=-1)
        else:
            alpha_sc = alpha_prev[..., 0].copy()
            # Undo the delta shift for n_s=1 case
            for n, (_, coeff) in enumerate(h_int):
                alpha_sc[n, 0, 0] -= -np.sign(coeff) * delta[0]
                alpha_sc[n, 0, 1] -= delta[1] * (abs(alpha_sc[n, 0, 1] - delta[1]) > 1e-6)
                alpha_sc[n, 1, 0] -= np.sign(coeff) * delta[1] * (abs(alpha_sc[n, 1, 0] - delta[1]) > 1e-6)
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
        gf_struct = self.constr_params['gf_struct']

        # Determine number of D0 alpha entries
        n_D0_total = 0
        if self.constr_params['use_D']:
            n_bl = len(gf_struct)
            R = gf_struct[0][1]
            n_D0_total = n_bl * n_bl * R * R

        assert solve_params['n_s'] == 2
        alpha = np.zeros((n_terms + n_D0_total, 2, 2, 2), dtype=complex if gtau_is_complex else float)
        for l, (term, _) in enumerate(h_int):
            bl0, bl1, u0, u0p, u1, u1p = self._indices_from_quartic_term(term)
            same_block = bl0 == bl1
            diag_u = u0 == u0p and u1 == u1p

            if not diag_u:
                if not same_block and u0 == u1p and u1 == u0p:
                    raise NotImplementedError("Spin-flip terms are not yet treated")
                if not same_block and u0 == u1 and u0p == u1p:
                    raise NotImplementedError("Pair-hopping terms are not yet treated")
                raise ValueError("Unknown term type")

            # density-density terms with diagonal u indices
            use_offdiag = same_block and u0 != u1
            d0, d1 = delta[0], delta[1] if use_offdiag else 0.0
            alpha[l, ..., 0] = [[0.5 + d0,  d1], [-d1, 0.5 - d0]]
            alpha[l, ..., 1] = [[0.5 - d0, -d1], [ d1, 0.5 + d0]]

        # Fill D0 alpha entries with trivial density (0.5) + delta shifts
        if n_D0_total > 0:
            for ibl1 in range(n_bl):
                for ibl2 in range(n_bl):
                    for a in range(R):
                        for b in range(R):
                            d = ibl1 * n_bl * R * R + ibl2 * R * R + a * R + b
                            alpha[n_terms + d, 0, 0, 0] = 0.5 + delta[0]
                            alpha[n_terms + d, 1, 1, 0] = 0.5 - delta[0]
                            alpha[n_terms + d, 0, 0, 1] = 0.5 - delta[0]
                            alpha[n_terms + d, 1, 1, 1] = 0.5 + delta[0]

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

            # Normalize delta to a 2-element list
            delta = solve_params.get('delta', [0.1, 0.1])
            if np.isscalar(delta):
                solve_params['delta'] = [delta, delta]
            elif len(delta) != 2:
                raise ValueError("delta must have exactly two components")

            alpha_mode = solve_params.pop('alpha_mode', "automatic")
            if alpha_mode == "automatic":
                alpha = self.find_alpha_from_HF_solver(solve_params)
            elif alpha_mode == "trivial":
                alpha = self.trivial_alpha(solve_params)
            else:
                raise ValueError(f"No such alpha_mode: {alpha_mode}")

            mpi_print(" --- Alpha Tensor : ")
            for s in range(solve_params['n_s']):
                if solve_params['n_s'] > 1:
                    mpi_print(f"Alpha Tensor s = {s + 1}:")
                mpi_print(str(alpha[..., s]))
            solve_params['alpha'] = alpha

        solve_status = SolverCore.solve(self, **solve_params)
        return solve_status
