# Copyright (c) 2018--present, The Simons Foundation
# This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

r"""User-facing CTINT solver.

This module exposes :class:`Solver`, a Python wrapper around
:class:`~triqs_ctint.solver_core.SolverCore`.
"""

from .solver_core import SolverCore, ConstrParamsT, SolveParamsT
from .version import gtau_is_complex, interaction_is_complex
from triqs.gfs import *
from triqs.utility import mpi
from triqs_hartree_fock import ImpuritySolver as HFSolver

import numpy as np


# === Alpha tensor validation and clipping

def _validate_and_clip_alpha(alpha):
    """
    Validate symmetry properties of the alpha tensor and clip values to the
    physical window. Each 2x2 alpha matrix corresponds to a 1-particle density
    matrix rho. Both rho and (1 - rho) must be positive definite, which requires
    diagonal elements in [0, 1] and |off-diagonal| <= 0.5.
    """
    eps = 1e-5
    clamp = lambda x, lo, hi: min(max(x, lo), hi)
    is_complex = np.iscomplexobj(alpha)

    for n in range(alpha.shape[0]):
        for s in range(alpha.shape[3]):
            a = alpha[n, :, :, s]

            if is_complex:
                if a[1, 0] != np.conj(a[0, 1]):
                    raise RuntimeError(f"Alpha tensor is not hermitian for term {n}, s={s}: "
                                       f"alpha[1,0]={a[1,0]}, conj(alpha[0,1])={np.conj(a[0, 1])}")
                for i in range(2):
                    if abs(a[i, i].imag) > 1e-12:
                        raise RuntimeError(f"Alpha tensor has non-real diagonal for term {n}, s={s}: "
                                           f"alpha[{i},{i}]={a[i, i]}")
                if abs(a[0, 1]) > 0.5 + eps:
                    raise RuntimeError(f"Alpha tensor |off-diagonal|={abs(a[0, 1]):.6f} > {0.5 + eps} "
                                       f"for term {n}, s={s}")
            else:
                if a[1, 0] != a[0, 1]:
                    raise RuntimeError(f"Alpha tensor is not symmetric for term {n}, s={s}: "
                                       f"alpha[1,0]={a[1, 0]}, alpha[0,1]={a[0, 1]}")
                a[0, 1] = clamp(a[0, 1], -0.5 - eps, 0.5 + eps)
                a[1, 0] = clamp(a[1, 0], -0.5 - eps, 0.5 + eps)

            a[0, 0] = clamp(a[0, 0].real, -eps, 1 + eps)
            a[1, 1] = clamp(a[1, 1].real, -eps, 1 + eps)

    return alpha


# === Some utility functions

def mpi_print(arg):
    """Print an object on the MPI master node only.

    Output from non-master ranks is suppressed so that a parallel run emits a
    single copy of the message. NumPy floating-point output is formatted with a
    precision of 4 digits for the duration of the call.

    Parameters
    ----------
    arg : object
        The object to print. It is passed unchanged to the built-in
        :func:`print`.

    Returns
    -------
    None
    """
    if mpi.is_master_node():
        with np.printoptions(precision=4):
            print(arg)

# === The SolverCore Wrapper

class Solver(SolverCore):
    r"""Continuous-time interaction-expansion impurity solver.

    Python wrapper around :class:`~triqs_ctint.solver_core.SolverCore`. Unlike a
    bare ``SolverCore``, this class can construct the ``alpha`` tensor
    automatically from a self-consistent Hartree-Fock solution when it is not
    given explicitly to :meth:`solve`.

    Parameters
    ----------
    **constr_params
        Construction parameters forwarded to
        :class:`~triqs_ctint.solver_core.ConstrParamsT`; see that class for the
        full list with defaults.
    """

    def __init__(self, **constr_params):
        """Initialise the solver. See :class:`Solver` for the constructor parameters."""
        constr_params['gf_struct'] = fix_gf_struct_type(constr_params['gf_struct'])

        # Initialise the core solver
        SolverCore.__init__(self, ConstrParamsT(**constr_params))


    def _indices_from_quartic_term(self, term):
        """Return (b0, b1, u0, u0p, u1, u1p) for a quartic term cdag cdag c c.

        :meta private:
        """
        bl0, u0 = term[0][1]
        bl1, u1 = term[1][1]
        bl1p, u1p = term[2][1]
        bl0p, u0p = term[3][1]
        assert bl0 == bl0p and bl1 == bl1p
        return (bl0, bl1, u0, u0p, u1, u1p)

    def find_alpha_from_HF_solver(self, solve_params):
        r"""Determine the :math:`\alpha`-tensor from a self-consistent Hartree-Fock solution.

        Runs :class:`triqs_hartree_fock.ImpuritySolver` on the same DLR mesh as
        the CT-INT solver to obtain the self-consistent Green's function
        :math:`G(i\omega)` and its density matrix :math:`\rho`. The
        :math:`\alpha`-tensor is built from the matrix elements of :math:`\rho`,
        with a :math:`\delta`-shift applied per auxiliary-spin component. If a
        previous solve is available, its :math:`\alpha`-tensor warm-starts the
        Hartree-Fock self-energy.

        Parameters
        ----------
        solve_params : dict
            The solve parameters. Must contain ``h_int``, the interaction
            Hamiltonian :math:`\hat H_\mathrm{int}`. The entries ``delta``
            (two-component shift, default ``[0.1, 0.1]``) and ``n_s`` (number of
            auxiliary spins, ``1`` or ``2``; default ``2``) are read; ``n_s`` is
            written back into ``solve_params``.

        Returns
        -------
        numpy.ndarray
            The :math:`\alpha`-tensor of shape
            ``(n_terms + n_D0, 2, 2, n_s)``, broadcast to all MPI ranks.
        """
        mpi_print("Determine alpha-tensor using triqs_hartree_fock")

        gf_struct = self.constr_params.gf_struct
        h_int = solve_params['h_int']
        delta = solve_params.pop('delta', [0.1, 0.1])
        n_s = solve_params.get('n_s', 2)
        assert n_s in [1, 2], "Solve parameter n_s has to be either 1 or 2 for automatic alpha mode"

        # The number of terms in h_int determines the leading dimension of alpha
        n_terms = len(list(h_int))

        # Create HF solver instance on the same DLR mesh as ctint
        hf_solver = HFSolver(
            gf_struct=gf_struct,
            mesh=self.G0_iw.mesh,
            dc=False,
            force_real=not gtau_is_complex
        )

        # Copy G0_iw to the HF solver (same DLR mesh)
        hf_solver.G0_iw << self.G0_iw

        # Initialize Sigma_HF from previous alpha if available
        if self.last_solve_params is not None:
            mpi_print("Initializing HF solver from previous iteration")
            alpha_prev = self.last_solve_params.alpha
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
        if self.constr_params.use_D:
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
                alpha[n, 1, 0, s] = np.conj(n01) + np.sign(coeff) * sgn * delta[1] * has_offdiag
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

        alpha = _validate_and_clip_alpha(alpha)
        alpha = mpi.bcast(alpha, root=0)

        # Make sure to set n_s as provided by the user
        solve_params['n_s'] = n_s

        return alpha

    def _initialize_hf_sigma_from_alpha(self, hf_solver, h_int, alpha_prev, delta):
        """Warm-start hf_solver.Sigma_HF from a previous alpha tensor (modified in place).

        :meta private:
        """
        gf_struct = self.constr_params.gf_struct

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
        r"""Build a simple :math:`\alpha`-tensor from a half-filling ansatz.

        Constructs the :math:`\alpha`-tensor for density-density interactions
        from a trivial density of :math:`1/2` per spin, shifted by
        :math:`\delta`, without solving any auxiliary problem. Requires
        ``n_s == 2``.

        Parameters
        ----------
        solve_params : dict
            The solve parameters (see :class:`~triqs_ctint.solver_core.SolveParamsT`). 
            Must contain ``h_int`` (the interaction Hamiltonian 
            :math:`\hat H_\mathrm{int}`) and ``n_s``, which must be ``2``. The entry 
            ``delta`` (two-component shift, default ``[0.5 + 1e-2, 1e-2]``) is read.

        Returns
        -------
        numpy.ndarray
            The :math:`\alpha`-tensor of shape ``(n_terms + n_D0, 2, 2, 2)``.

        Raises
        ------
        NotImplementedError
            If the interaction contains spin-flip or pair-hopping terms.
        ValueError
            If a term has an unrecognised structure.
        """
        h_int = solve_params['h_int']
        n_terms = len(list(h_int))
        delta = solve_params.pop('delta', [0.5 + 1e-2, 1e-2])
        gf_struct = self.constr_params.gf_struct

        # Determine number of D0 alpha entries
        n_D0_total = 0
        if self.constr_params.use_D:
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

        return _validate_and_clip_alpha(alpha)

    def solve(self, **solve_params):
        r"""
        Solve the impurity problem.

        Parameters
        ----------
        **solve_params
            Solve parameters forwarded to :class:`~triqs_ctint.solver_core.SolveParamsT`,
            which is passed to :meth:`~triqs_ctint.solver_core.SolverCore.solve`.
            The only two required parameters are (i) ``h_int``, the local interaction 
            Hamiltonian :math:`\hat H_\mathrm{int}`, and (ii) ``n_cycles``, the number 
            of Monte-Carlo cycles. Note that the :math:`\alpha` tensor is optional. If 
            it is not given, it is constructed from the density matrix of the 
            self-consistent Hartree-Fock solution.
        delta : float, optional
            Value of :math:`\delta` used to construct the :math:`\alpha`-tensor
            (default ``0.1``). A larger :math:`\delta` improves the Monte-Carlo
            sign at the cost of a larger perturbation order. Only used if
            ``alpha`` is not given explicitly.

        Returns
        -------
        solve_status
            The Monte-Carlo run status returned by the core
            :meth:`~triqs_ctint.solver_core.SolverCore.solve`.
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

        solve_status = SolverCore.solve(self, SolveParamsT(**solve_params))
        return solve_status
