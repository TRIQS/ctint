"""Verification of chirality-chirality correlators <T*T> against exact diagonalization.

Measures the real-valued chirality-chirality correlators T_ijk * T_lmn (degree 12)
on a 2x2 Hubbard plaquette. The individual <T_ijk> vanish by complex-conjugation
symmetry of the real Hamiltonian, but <T*T> is nonzero.

Exercises compute_insertk_ratio for k up to 6 per block.

Site layout:
    0 - 1
    |   |
    2 - 3
"""

from triqs_ctint import Solver
from triqs.gf import *
from triqs.operators import n, c_dag, c
from triqs.atom_diag import AtomDiagComplex, atomic_density_matrix, trace_rho_op
import numpy as np
import triqs.utility.mpi as mpi

# === Physical parameters: 2x2 plaquette with broken symmetry ===
beta = 5.0
N_sites = 4

eps = [-0.8, -1.1, -0.5, -0.7]        # site-dependent on-site energy
h_z = [0.15, -0.10, 0.08, -0.05]      # Zeeman field (breaks time-reversal)
U_vals = [1.5, 2.0, 1.2, 1.8]         # site-dependent Hubbard U
t_01, t_23, t_02, t_13 = 0.5, 0.3, 0.4, 0.35  # asymmetric hoppings

gf_struct = [("up", N_sites), ("dn", N_sites)]

# Build hopping matrix (same for up/dn before Zeeman)
t_mat = np.zeros((N_sites, N_sites))
t_mat[0, 1] = t_mat[1, 0] = -t_01
t_mat[2, 3] = t_mat[3, 2] = -t_23
t_mat[0, 2] = t_mat[2, 0] = -t_02
t_mat[1, 3] = t_mat[3, 1] = -t_13

hloc0_up = t_mat.copy()
hloc0_dn = t_mat.copy()
for i in range(N_sites):
    hloc0_up[i, i] = eps[i] + h_z[i]
    hloc0_dn[i, i] = eps[i] - h_z[i]

h_int = sum(U_vals[i] * n("up", i) * n("dn", i) for i in range(N_sites))


# === Spin operators ===
def S_x(i):
    return 0.5 * (c_dag("up", i) * c("dn", i) + c_dag("dn", i) * c("up", i))


def S_y(i):
    return -0.5j * (c_dag("up", i) * c("dn", i) - c_dag("dn", i) * c("up", i))


def S_z(i):
    return 0.5 * (n("up", i) - n("dn", i))


def chirality(i, j, k):
    """Scalar spin chirality T_ijk = S_i . (S_j x S_k)"""
    return (
        S_x(i) * (S_y(j) * S_z(k) - S_z(j) * S_y(k))
        + S_y(i) * (S_z(j) * S_x(k) - S_x(j) * S_z(k))
        + S_z(i) * (S_x(j) * S_y(k) - S_y(j) * S_x(k))
    )


# === Exact diagonalization reference ===
fops = [("dn", i) for i in range(N_sites)] + [("up", i) for i in range(N_sites)]
H_full = sum(eps[i] * (n("up", i) + n("dn", i)) + h_z[i] * (n("up", i) - n("dn", i)) for i in range(N_sites))
H_full += sum(U_vals[i] * n("up", i) * n("dn", i) for i in range(N_sites))
for s in ["up", "dn"]:
    H_full += -t_01 * (c_dag(s, 0) * c(s, 1) + c_dag(s, 1) * c(s, 0))
    H_full += -t_23 * (c_dag(s, 2) * c(s, 3) + c_dag(s, 3) * c(s, 2))
    H_full += -t_02 * (c_dag(s, 0) * c(s, 2) + c_dag(s, 2) * c(s, 0))
    H_full += -t_13 * (c_dag(s, 1) * c(s, 3) + c_dag(s, 3) * c(s, 1))

ad = AtomDiagComplex(H_full, fops)
dm = atomic_density_matrix(ad, beta)


def ed_expect(op):
    return trace_rho_op(dm, op, ad)


# === Define test operators ===
T_012 = chirality(0, 1, 2)
T_013 = chirality(0, 1, 3)
T_123 = chirality(1, 2, 3)

test_ops = {
    # Degree 12: chirality-chirality correlators (k up to 6 per block, real-valued)
    "T_012*T_012": T_012 * T_012,
    "T_012*T_123": T_012 * T_123,
    "T_012*T_013": T_012 * T_013,
}

ed_values = {name: complex(ed_expect(op)).real for name, op in test_ops.items()}

# === chiAB operator pairs ===
chi_pairs = {
    "chiAB_T012_T012": (T_012, T_012),
}

chi_expected = {
    "chiAB_T012_T012": ed_values["T_012*T_012"],
}

chi_to_static = {
    "chiAB_T012_T012": 0,  # index of T_012*T_012 in test_ops
}

# === Run solver ===
S = Solver(beta=beta, gf_struct=gf_struct, n_tau=201, dlr_wmax=5.0, dlr_eps=1e-7)

S.G0_iw["up"] << inverse(iOmega_n - hloc0_up)
S.G0_iw["dn"] << inverse(iOmega_n - hloc0_dn)

S.solve(
    h_int=h_int,
    delta=0.5,
    n_cycles=20000,
    length_cycle=50,
    n_warmup_cycles=1000,
    measure_static_obs=True,
    static_obs=list(test_ops.values()),
    n_tau_static_obs=10,
    measure_chiAB_tau=True,
    chi_ops=list(chi_pairs.values()),
    random_seed=54321,
)

# === Verification ===
if mpi.is_master_node():
    tol = 0.03
    all_passed = True
    op_names = list(test_ops.keys())
    chi_names = list(chi_pairs.keys())

    print("=" * 70)
    print("STATIC_OBS vs EXACT DIAGONALIZATION")
    print("=" * 70)
    for i, name in enumerate(op_names):
        mc_val = S.static_obs[i].real
        ed_val = ed_values[name]
        mc_err = S.static_obs_errors[i]
        diff = mc_val - ed_val
        ok = abs(diff) < max(tol, 3 * mc_err)
        sign_ok = (mc_val > 0) == (ed_val > 0) or abs(ed_val) < tol
        status = "PASS" if (ok and sign_ok) else "FAIL"
        if not (ok and sign_ok):
            all_passed = False
        print(f"  {name:20s}  MC={mc_val:+.6f}  ED={ed_val:+.6f}  diff={diff:+.2e}  err={mc_err:.2e}  [{status}]")

    print()
    print("=" * 70)
    print("CHIAB_TAU(tau=0) vs EXACT DIAGONALIZATION")
    print("=" * 70)
    for i, name in enumerate(chi_names):
        chi_tau0 = S.chiAB_tau.data[0, i].real
        ed_val = chi_expected[name]
        diff = chi_tau0 - ed_val
        ok = abs(diff) < tol
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_passed = False
        print(f"  {name:20s}  MC={chi_tau0:+.6f}  ED={ed_val:+.6f}  diff={diff:+.2e}  [{status}]")

    print()
    print("=" * 70)
    print("CHIAB_TAU(tau=0) vs STATIC_OBS (cross-validation)")
    print("=" * 70)
    for name, static_idx in chi_to_static.items():
        i = chi_names.index(name)
        chi_val = S.chiAB_tau.data[0, i].real
        static_val = S.static_obs[static_idx].real
        diff = chi_val - static_val
        ok = abs(diff) < tol
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_passed = False
        print(f"  {name:20s}  chiAB={chi_val:+.6f}  static={static_val:+.6f}  diff={diff:+.2e}  [{status}]")

    print()
    if all_passed:
        print("All verification tests PASSED!")
    else:
        raise AssertionError("Some verification tests FAILED -- see details above")
