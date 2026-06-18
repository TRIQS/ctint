"""Test static_obs measurement against densities and chiAB_tau."""

from triqs_ctint import Solver
from triqs.gf import *
from triqs.operators import n, c_dag, c
import numpy as np
import triqs.utility.mpi as mpi

test_name = "static_obs"

# Physical parameters (half-filled 2-site Anderson model)
beta = 10.0
U = 2.0

gf_struct = [("up", 2), ("dn", 2)]

# Hopping between sites 0 and 1
t = 0.5
hloc0 = np.array([[-U / 2, -t], [-t, -U / 2]])

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_tau=201, dlr_wmax=10.0, dlr_eps=1e-10)

# Initialize G0
for bl, g_bl in S.G0_iw:
    g_bl << inverse(iOmega_n - hloc0)

# Interaction
h_int = U * n("up", 0) * n("dn", 0) + U * n("up", 1) * n("dn", 1)

# Define static observables
static_ops = [
    n("up", 0),                          # density (bilinear)
    n("dn", 0),                          # density (bilinear)
    n("up", 0) * n("dn", 0),            # double occupancy (quartic, AAAA)
    n("up", 0) * n("dn", 1),            # inter-block density-density (quartic, AABB)
    n("up", 0) + n("dn", 0),            # sum of bilinears
]

# Also measure chiAB_tau for cross-validation at tau=0
chi_ops = [
    (n("up", 0), n("dn", 0)),           # same as static_obs[2]
    (n("up", 0), n("dn", 1)),           # same as static_obs[3]
]

S.solve(h_int=h_int,
        delta=0.5,
        n_cycles=5000,
        length_cycle=50,
        n_warmup_cycles=500,
        measure_densities=True,
        measure_chiAB_tau=True,
        chi_ops=chi_ops,
        measure_static_obs=True,
        static_obs=static_ops,
        n_tau_static_obs=10)

if mpi.is_master_node():
    # Test 1: Bilinear vs densities
    n_up_0_static = S.static_obs[0].real
    n_dn_0_static = S.static_obs[1].real
    # The solver preserves the given gf_struct block order [("up", 2), ("dn", 2)]
    n_up_0_dens = S.densities[0][0]  # "up" block is index 0
    n_dn_0_dens = S.densities[1][0]  # "dn" block is index 1

    print(f"n_up_0 from static_obs: {n_up_0_static:.10f}")
    print(f"n_up_0 from densities:  {n_up_0_dens:.10f}")
    print(f"n_up_0 diff:            {n_up_0_static - n_up_0_dens:.2e}")
    print(f"n_dn_0 from static_obs: {n_dn_0_static:.10f}")
    print(f"n_dn_0 from densities:  {n_dn_0_dens:.10f}")
    print(f"n_dn_0 diff:            {n_dn_0_static - n_dn_0_dens:.2e}")
    assert abs(n_up_0_static - n_up_0_dens) < 0.05, f"Bilinear mismatch n_up: {n_up_0_static} vs {n_up_0_dens}"
    assert abs(n_dn_0_static - n_dn_0_dens) < 0.05, f"Bilinear mismatch n_dn: {n_dn_0_static} vs {n_dn_0_dens}"

    # Test 2: Sum of bilinears
    n_tot_0 = S.static_obs[4].real
    n_up_plus_dn = S.static_obs[0].real + S.static_obs[1].real
    print(f"n_up_0 + n_dn_0 from static_obs: {n_tot_0:.6f}")
    print(f"sum of individual:                {n_up_plus_dn:.6f}")
    assert abs(n_tot_0 - n_up_plus_dn) < 0.05, f"Sum mismatch: {n_tot_0} vs {n_up_plus_dn}"

    # Test 3: Quartic (AAAA) vs chiAB_tau at tau=0
    dd_00_static = S.static_obs[2].real
    # chiAB_tau data at tau=0, operator index 0
    dd_00_chi = S.chiAB_tau.data[0, 0].real
    print(f"n_up_0*n_dn_0 from static_obs: {dd_00_static:.6f}")
    print(f"n_up_0*n_dn_0 from chiAB_tau:  {dd_00_chi:.6f}")
    assert abs(dd_00_static - dd_00_chi) < 0.05, f"Quartic AAAA mismatch: {dd_00_static} vs {dd_00_chi}"

    # Test 4: Quartic (AABB) vs chiAB_tau at tau=0
    dd_01_static = S.static_obs[3].real
    dd_01_chi = S.chiAB_tau.data[0, 1].real
    print(f"n_up_0*n_dn_1 from static_obs: {dd_01_static:.6f}")
    print(f"n_up_0*n_dn_1 from chiAB_tau:  {dd_01_chi:.6f}")
    assert abs(dd_01_static - dd_01_chi) < 0.05, f"Quartic AABB mismatch: {dd_01_static} vs {dd_01_chi}"

    print("\nAll static_obs tests passed!")
