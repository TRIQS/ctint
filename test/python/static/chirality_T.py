"""Verification of scalar spin chirality <T_ijk> against exact diagonalization.

Measures the scalar spin chirality T_ijk = S_i . (S_j x S_k) (degree 6) on a
2x2 Hubbard plaquette with spin-dependent magnetic flux (Peierls phases).
Up-spins see flux +phi, dn-spins see -phi around the plaquette. This mimics
spin-orbit-induced effective flux and breaks time-reversal symmetry while
generating non-coplanar spin correlations that make <T_ijk> nonzero (~0.02).

A spin-independent flux alone gives negligible chirality because both spin
species carry the same orbital current. Spin-dependent flux creates differential
spin currents that directly couple to the scalar chirality.

Requires GTAU_IS_COMPLEX=ON at compile time for complex G(tau) support.

Uses a single-block gf_struct to accommodate the spin-dependent complex hoppings.
Orbital layout within the block: 0..3 = up (sites 0-3), 4..7 = dn (sites 0-3).

Site layout:
    0 - 1
    |   |
    2 - 3
"""

from triqs_ctint import Solver, version
from triqs.gf import *
from triqs.operators import n, c_dag, c
from triqs.atom_diag import AtomDiagComplex, atomic_density_matrix, trace_rho_op
import numpy as np
import triqs.utility.mpi as mpi

if not version.gtau_is_complex:
    print("FIXME: GTAU_IS_COMPLEX not enabled, skipping test")
    exit()

# === Orbital mapping: single block "s" with 8 orbitals ===
# Orbitals 0..3 = spin up at sites 0..3
# Orbitals 4..7 = spin dn at sites 0..3
N_sites = 4
N_orb = 2 * N_sites
BL = "s"


def up(site):
    return site


def dn(site):
    return site + N_sites


# === Physical parameters: 2x2 plaquette with spin-dependent flux ===
beta = 5.0

eps = [-0.8, -1.1, -0.5, -0.7]        # site-dependent on-site energy
h_z = [0.15, -0.10, 0.08, -0.05]      # Zeeman field
U_vals = [1.5, 2.0, 1.2, 1.8]         # site-dependent Hubbard U

# Hopping magnitudes around the plaquette (0→1→3→2→0)
t_mag = {(0, 1): 0.5, (1, 3): 0.35, (3, 2): 0.3, (2, 0): 0.4}

# Spin-dependent Peierls phase: up-spins see +phi, dn-spins see -phi.
# Phase distributed uniformly around the plaquette: each bond gets phi/4.
phi = np.pi  # total flux (pi gives max chirality ~0.02)

gf_struct = [(BL, N_orb)]

# Build 8x8 hopping matrix with spin-dependent Peierls phases
hloc0 = np.zeros((N_orb, N_orb), dtype=complex)

# On-site: eps + Zeeman
for i in range(N_sites):
    hloc0[up(i), up(i)] = eps[i] + h_z[i]
    hloc0[dn(i), dn(i)] = eps[i] - h_z[i]

# Spin-dependent hopping: up gets phase +phi/4, dn gets -phi/4 per bond
p = phi / 4
for (i, j), t in t_mag.items():
    t_up = t * np.exp(1j * p)
    hloc0[up(i), up(j)] += -t_up
    hloc0[up(j), up(i)] += -np.conj(t_up)
    t_dn = t * np.exp(-1j * p)
    hloc0[dn(i), dn(j)] += -t_dn
    hloc0[dn(j), dn(i)] += -np.conj(t_dn)

# Verify Hermiticity
assert np.allclose(hloc0, hloc0.conj().T), "hloc0 is not Hermitian"

# Interaction
h_int = sum(U_vals[i] * n(BL, up(i)) * n(BL, dn(i)) for i in range(N_sites))


# === Spin operators (single-block) ===
def S_x(i):
    return 0.5 * (c_dag(BL, up(i)) * c(BL, dn(i)) + c_dag(BL, dn(i)) * c(BL, up(i)))


def S_y(i):
    return -0.5j * (c_dag(BL, up(i)) * c(BL, dn(i)) - c_dag(BL, dn(i)) * c(BL, up(i)))


def S_z(i):
    return 0.5 * (n(BL, up(i)) - n(BL, dn(i)))


def chirality(i, j, k):
    """Scalar spin chirality T_ijk = S_i . (S_j x S_k)"""
    return (
        S_x(i) * (S_y(j) * S_z(k) - S_z(j) * S_y(k))
        + S_y(i) * (S_z(j) * S_x(k) - S_x(j) * S_z(k))
        + S_z(i) * (S_x(j) * S_y(k) - S_y(j) * S_x(k))
    )


# === Exact diagonalization reference ===
# Build H_full using the same single-block operators
fops = [(BL, i) for i in range(N_orb)]

H_full = sum(eps[i] * (n(BL, up(i)) + n(BL, dn(i))) + h_z[i] * (n(BL, up(i)) - n(BL, dn(i))) for i in range(N_sites))
H_full += sum(U_vals[i] * n(BL, up(i)) * n(BL, dn(i)) for i in range(N_sites))

# Spin-dependent hoppings with Peierls phases
for (i, j), t in t_mag.items():
    t_up = t * np.exp(1j * p)
    H_full += -t_up * c_dag(BL, up(i)) * c(BL, up(j)) - np.conj(t_up) * c_dag(BL, up(j)) * c(BL, up(i))
    t_dn = t * np.exp(-1j * p)
    H_full += -t_dn * c_dag(BL, dn(i)) * c(BL, dn(j)) - np.conj(t_dn) * c_dag(BL, dn(j)) * c(BL, dn(i))

ad = AtomDiagComplex(H_full, fops)
dm = atomic_density_matrix(ad, beta)


def ed_expect(op):
    return trace_rho_op(dm, op, ad)


# === Define test operators ===
T_012 = chirality(0, 1, 2)
T_013 = chirality(0, 1, 3)
T_023 = chirality(0, 2, 3)
T_123 = chirality(1, 2, 3)

test_ops = {
    "T_012": T_012,
    "T_013": T_013,
    "T_023": T_023,
    "T_123": T_123,
}

# <T_ijk> is generally complex with SOC; for ED comparison use full complex value
ed_values = {name: complex(ed_expect(op)) for name, op in test_ops.items()}

# === chiAB operator pairs: T(tau) * T(0) at tau=0 gives <T*T> ===
chi_pairs = {
    "chiAB_T012_T012": (T_012, T_012),
    "chiAB_T012_T123": (T_012, T_123),
}

chi_products = {
    "chiAB_T012_T012": T_012 * T_012,
    "chiAB_T012_T123": T_012 * T_123,
}
chi_expected = {name: complex(ed_expect(op)) for name, op in chi_products.items()}

# === Run solver ===
S = Solver(beta=beta, gf_struct=gf_struct, n_tau=201, dlr_wmax=5.0, dlr_eps=1e-7)

S.G0_iw[BL] << inverse(iOmega_n - hloc0)

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
    print("STATIC_OBS <T_ijk> vs EXACT DIAGONALIZATION")
    print("=" * 70)
    for i, name in enumerate(op_names):
        mc_val = complex(S.static_obs[i])
        ed_val = ed_values[name]
        mc_err = S.static_obs_errors[i]
        diff = abs(mc_val - ed_val)
        ok = diff < max(tol, 3 * mc_err)
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_passed = False
        print(f"  {name:10s}  MC={mc_val:+.6f}  ED={ed_val:+.6f}  |diff|={diff:.2e}  err={mc_err:.2e}  [{status}]")

    print()
    print("=" * 70)
    print("CHIAB_TAU(tau=0) <T*T> vs EXACT DIAGONALIZATION")
    print("=" * 70)
    for i, name in enumerate(chi_names):
        chi_tau0 = complex(S.chiAB_tau.data[0, i])
        ed_val = chi_expected[name]
        diff = abs(chi_tau0 - ed_val)
        ok = diff < tol
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_passed = False
        print(f"  {name:25s}  MC={chi_tau0:+.6f}  ED={ed_val:+.6f}  |diff|={diff:.2e}  [{status}]")

    print()
    if all_passed:
        print("All verification tests PASSED!")
    else:
        raise AssertionError("Some verification tests FAILED -- see details above")
