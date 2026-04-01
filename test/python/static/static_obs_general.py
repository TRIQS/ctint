"""Comprehensive verification of arbitrary-degree operator measurements.

Compares static_obs and chiAB_tau results against exact diagonalization (ED)
for a 3-site Anderson model with broken site and spin symmetry, covering
all code paths:
- Bilinear (degree 2, k=1 per block)
- Quartic AABB (degree 4, k=1+k=1 across blocks)
- Quartic AAAA (degree 4, k=2 in one block)
- Sextic mixed (degree 6, k=2+k=1)
- Sextic AAAA (degree 6, k=3 in one block) -- exercises compute_insertk_ratio
- chiAB cross-validation at tau~0 for AABB, AAAA, and sextic pairs
- chiAB with unbalanced A/B monomials (A has more c-dags than c's, B vice versa)
"""

from triqs_ctint import Solver
from triqs.gf import *
from triqs.operators import n, c_dag, c
from triqs.atom_diag import AtomDiag, atomic_density_matrix, trace_rho_op
import numpy as np
import triqs.utility.mpi as mpi

# === Physical parameters: 3-site model with broken site + spin symmetry ===
#
# H = Σ_i [ ε_i (n_up,i + n_dn,i) + h_i (n_up,i - n_dn,i) + U_i n_up,i n_dn,i ]
#     - t_01 Σ_σ (c†_σ,0 c_σ,1 + h.c.)
#     - t_12 Σ_σ (c†_σ,1 c_σ,2 + h.c.)
#
beta = 10.0
eps = [-0.8, -1.1, -0.5]   # site-dependent local energy
h_z = [0.05, -0.10, 0.03]  # site-dependent Zeeman field
U_vals = [2.0, 3.0, 1.5]   # site-dependent Hubbard U
t01, t12 = 0.5, 0.3        # asymmetric hoppings (no t02)

gf_struct = [("up", 3), ("dn", 3)]

# h_loc for G0: on-site energy = eps_i + sigma * h_i, plus hopping
hloc0_up = np.array([[eps[0] + h_z[0], -t01, 0], [-t01, eps[1] + h_z[1], -t12], [0, -t12, eps[2] + h_z[2]]])
hloc0_dn = np.array([[eps[0] - h_z[0], -t01, 0], [-t01, eps[1] - h_z[1], -t12], [0, -t12, eps[2] - h_z[2]]])

h_int = sum(U_vals[i] * n("up", i) * n("dn", i) for i in range(3))

# === Exact diagonalization reference ===
fops = [("dn", i) for i in range(3)] + [("up", i) for i in range(3)]
H_full = sum(eps[i] * (n("up", i) + n("dn", i)) + h_z[i] * (n("up", i) - n("dn", i)) for i in range(3))
H_full += sum(U_vals[i] * n("up", i) * n("dn", i) for i in range(3))
H_full += -t01 * sum(c_dag(s, 0) * c(s, 1) + c_dag(s, 1) * c(s, 0) for s in ["up", "dn"])
H_full += -t12 * sum(c_dag(s, 1) * c(s, 2) + c_dag(s, 2) * c(s, 1) for s in ["up", "dn"])

ad = AtomDiag(H_full, fops)
dm = atomic_density_matrix(ad, beta)


def ed_expect(op):
    return trace_rho_op(dm, op, ad)


# === Define ALL test operators ===
test_ops = {
    # Bilinear (degree 2, k=1)
    "n_up(0)": n("up", 0),
    "n_dn(0)": n("dn", 0),
    "n_up(1)": n("up", 1),
    # Quartic AABB (degree 4, k=1+k=1 cross-block)
    "n_up(0)*n_dn(0)": n("up", 0) * n("dn", 0),
    "n_up(0)*n_dn(1)": n("up", 0) * n("dn", 1),
    # Quartic AAAA (degree 4, k=2 in one block -> insert2_ratios)
    "n_up(0)*n_up(1)": n("up", 0) * n("up", 1),
    "n_dn(0)*n_dn(1)": n("dn", 0) * n("dn", 1),
    # Sextic mixed (degree 6, k=2+k=1)
    "n_up(0)*n_dn(0)*n_up(1)": n("up", 0) * n("dn", 0) * n("up", 1),
    "n_up(0)*n_dn(0)*n_dn(1)": n("up", 0) * n("dn", 0) * n("dn", 1),
    # Sextic AAAA (degree 6, k=3 in one block -> compute_insertk_ratio)
    "n_up(0)*n_up(1)*n_up(2)": n("up", 0) * n("up", 1) * n("up", 2),
    "n_dn(0)*n_dn(1)*n_dn(2)": n("dn", 0) * n("dn", 1) * n("dn", 2),
    # Sextic product of unbalanced monomials (A*B balanced, but A and B individually are not)
    "cd_up0*cd_dn0*c_up1 * cd_up1*c_dn0*c_up0":
        c_dag("up", 0) * c_dag("dn", 0) * c("up", 1) * c_dag("up", 1) * c("dn", 0) * c("up", 0),
}

# Compute ED references
ed_values = {name: ed_expect(op).real for name, op in test_ops.items()}

# === chiAB operator pairs ===
chi_pairs = {
    # Bilinear x Bilinear -> quartic AABB
    "chiAB_AABB": (n("up", 0), n("dn", 0)),
    # Bilinear x Bilinear -> quartic AAAA (both in same block)
    "chiAB_AAAA": (n("up", 0), n("up", 1)),
    # Quartic x Bilinear -> sextic mixed
    "chiAB_sextic_mixed": (n("up", 0) * n("dn", 0), n("up", 1)),
    # Bilinear x Bilinear -> quartic AAAA for k=3 cross-check:
    # chiAB with (n_up(0)*n_up(1), n_up(2)) tests sextic AAAA in chiAB path
    "chiAB_sextic_k3": (n("up", 0) * n("up", 1), n("up", 2)),
    # Unbalanced monomials: A has 2 c-dags + 1 c, B has 1 c-dag + 2 c's.
    # Product A*B is balanced (3+3). Regression test for sign factor (-1)^{n_B_cdag * m_A_c}.
    "chiAB_unbalanced": (c_dag("up", 0) * c_dag("dn", 0) * c("up", 1),
                         c_dag("up", 1) * c("dn", 0) * c("up", 0)),
}

# chiAB at tau~0 should equal the product operator's expectation value
chi_expected = {
    "chiAB_AABB": ed_values["n_up(0)*n_dn(0)"],
    "chiAB_AAAA": ed_values["n_up(0)*n_up(1)"],
    "chiAB_sextic_mixed": ed_values["n_up(0)*n_dn(0)*n_up(1)"],
    "chiAB_sextic_k3": ed_values["n_up(0)*n_up(1)*n_up(2)"],
    "chiAB_unbalanced": ed_values["cd_up0*cd_dn0*c_up1 * cd_up1*c_dn0*c_up0"],
}

# === Run solver ===
S = Solver(beta=beta, gf_struct=gf_struct, n_tau=201, dlr_wmax=10.0, dlr_eps=1e-10)

S.G0_iw["up"] << inverse(iOmega_n - hloc0_up)
S.G0_iw["dn"] << inverse(iOmega_n - hloc0_dn)

S.solve(
    h_int=h_int,
    delta=0.1,
    n_cycles=80000,
    measure_static_obs=True,
    static_obs=list(test_ops.values()),
    n_tau_static_obs=10,
    measure_chiAB_tau=True,
    chi_ops=list(chi_pairs.values()),
    random_seed=98765,
)

# === Verification ===
if mpi.is_master_node():
    tol = 0.002
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
        print(f"  {name:30s}  MC={mc_val:+.6f}  ED={ed_val:+.6f}  diff={diff:+.2e}  err={mc_err:.2e}  [{status}]")

    print()
    print("=" * 70)
    print("CHIAB_TAU(tau~0) vs ED")
    print("=" * 70)
    for i, name in enumerate(chi_names):
        chi_tau0 = S.chiAB_tau.data[0, i].real
        ed_val = chi_expected[name]
        diff = chi_tau0 - ed_val
        ok = abs(diff) < tol
        sign_ok = (chi_tau0 > 0) == (ed_val > 0) or abs(ed_val) < tol
        status = "PASS" if (ok and sign_ok) else "FAIL"
        if not (ok and sign_ok):
            all_passed = False
        print(f"  {name:30s}  MC={chi_tau0:+.6f}  ED={ed_val:+.6f}  diff={diff:+.2e}  [{status}]")

    print()
    print("=" * 70)
    print("CHIAB_TAU(tau~0) vs STATIC_OBS (internal cross-validation)")
    print("=" * 70)
    # Map each chiAB pair to the corresponding static_obs index
    chi_to_static = {"chiAB_AABB": 3, "chiAB_AAAA": 5, "chiAB_sextic_mixed": 7, "chiAB_sextic_k3": 9, "chiAB_unbalanced": 11}
    for name, static_idx in chi_to_static.items():
        i = chi_names.index(name)
        chi_val = S.chiAB_tau.data[0, i].real
        static_val = S.static_obs[static_idx].real
        diff = chi_val - static_val
        ok = abs(diff) < tol
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_passed = False
        print(f"  {name:30s}  chiAB={chi_val:+.6f}  static={static_val:+.6f}  diff={diff:+.2e}  [{status}]")

    print()
    if all_passed:
        print("All verification tests PASSED!")
    else:
        raise AssertionError("Some verification tests FAILED -- see details above")
