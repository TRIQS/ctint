"""
Generic solver interface for CT-INT.

Provides a functional API to the triqs_ctint solver accepting
Delta_iw and h_loc0 (hybridization + one-body Hamiltonian) as input,
consistent with the cthyb and ctseg generic interfaces.
"""

import triqs.utility.mpi as mpi

from triqs.gfs import (
    MeshDLRImFreq, MeshImFreq,
    iOmega_n,
)
from triqs.gfs.tools import inverse
from triqs.solver_utils import SolverResults, make_gf_dlr_imfreq

from triqs_ctint import Solver


# Constructor parameter names for Solver (keyword-only)
_CONSTR_PARAM_NAMES = {
    'beta', 'gf_struct', 'dlr_wmax', 'dlr_eps', 'n_tau',
    'use_D', 'use_Jperp', 'n_tau_dynamical_interactions',
    'gtau_is_complex', 'interaction_is_complex',
}


def _build_G0_iw(Delta_iw, h_loc0):
    """Compute G0_iw = (iw - h_loc0 - Delta_iw)^{-1} block by block."""
    G0_iw = Delta_iw.copy()
    for ibl, (bl, gf) in enumerate(G0_iw):
        gf << inverse(iOmega_n - h_loc0[ibl] - Delta_iw[bl])
    return G0_iw


def solve_generic(
    Delta_iw,
    h_loc0,
    h_int,
    **solver_params,
):
    """Solve the impurity problem using CT-INT.

    Parameters
    ----------
    Delta_iw : BlockGf
        Hybridization function on MeshImFreq or MeshDLRImFreq.
    h_loc0 : list of np.ndarray
        One-body Hamiltonian as per-block matrices.
    h_int : Operator
        Interaction Hamiltonian.
    **solver_params
        Parameters forwarded to Solver constructor and solve().
        Constructor params: beta, gf_struct, dlr_wmax, dlr_eps, n_tau, use_D, use_Jperp,
                           gtau_is_complex, interaction_is_complex.
        Solve params: n_cycles, alpha, alpha_mode, delta, and all other SolverCore solve params.

    Returns
    -------
    SolverResults
    """
    mesh = Delta_iw.mesh
    params = solver_params.copy()

    # Extract gf_struct from Delta_iw if not provided
    if 'gf_struct' not in params:
        params['gf_struct'] = [(bl, gf.target_shape[0]) for (bl, gf) in Delta_iw]
    gf_struct = params['gf_struct']

    # Extract beta
    beta = params.pop('beta', mesh.beta)

    # Build G0_iw on the input mesh
    G0_iw = _build_G0_iw(Delta_iw, h_loc0)

    # If input is on MeshImFreq, project to DLR
    if isinstance(mesh, MeshImFreq):
        dlr_wmax = params.get('dlr_wmax')
        dlr_eps = params.get('dlr_eps', 1e-10)
        if dlr_wmax is None:
            raise ValueError("dlr_wmax must be provided when Delta_iw is on MeshImFreq")
        G0_iw, _ = make_gf_dlr_imfreq(G0_iw, w_max=dlr_wmax, eps=dlr_eps)
    elif isinstance(mesh, MeshDLRImFreq):
        # Extract DLR parameters from mesh if not provided
        params.setdefault('dlr_wmax', mesh.w_max)
        params.setdefault('dlr_eps', mesh.eps)
    else:
        raise NotImplementedError(f"Unsupported mesh type: {type(mesh)}")

    # Split into constructor and solve parameters
    constr_params = {'beta': beta, 'gf_struct': gf_struct}
    solve_kw = {}
    for key, val in params.items():
        if key in _CONSTR_PARAM_NAMES:
            constr_params[key] = val
        else:
            solve_kw[key] = val

    # Create solver and set G0_iw
    S = Solver(**constr_params)
    S.G0_iw << G0_iw

    # Solve
    mpi.report("Solving the impurity problem with CT-INT")
    solve_kw['h_int'] = h_int
    S.solve(**solve_kw)

    # Build results
    # Full Sigma = Sigma_dyn + Sigma_hartree
    Sigma_iw = S.Sigma_dyn_iw.copy()
    if S.Sigma_hartree is not None:
        for ibl, (bl, gf) in enumerate(Sigma_iw):
            gf += S.Sigma_hartree[ibl]

    result_kwargs = dict(
        G_iw=S.G_iw,
        Sigma_iw=Sigma_iw,
        Sigma_dynamic=S.Sigma_dyn_iw,
        Solver=S,
    )

    if S.Sigma_hartree is not None:
        result_kwargs['Sigma_HartreeFock'] = list(S.Sigma_hartree)

    if S.density_matrix is not None:
        result_kwargs['density_matrix'] = S.density_matrix

    return SolverResults(**result_kwargs)
