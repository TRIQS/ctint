"""
Test anomalous chiAB measurement (c†c† and cc operators) against exact atomic result.

Measures the s-wave pairing susceptibility chi_pair(tau) = <Delta†(tau) Delta(0)>
where Delta† = c†_up c†_dn (pair creation) and Delta = c_dn c_up (pair annihilation).

Uses a single-site attractive Hubbard model (atomic limit) where the exact result is known:
  chi_pair(tau) = exp(E_d * (tau - beta)) / Z
  with E_d = -2*mu + U, Z = 1 + 2*exp(-beta*E_1) + exp(-beta*E_d), E_1 = -mu
"""

from triqs_ctint import Solver
from triqs.gf import *
from triqs.operators import *
import triqs.utility.mpi as mpi
import numpy as np

# Physical parameters
U = -2.0
mu = U / 4.0  # Away from half-filling for nontrivial tau dependence
beta = 10.0

# MC parameters
n_cycles = 10000
length_cycle = 50
n_warmup = 500
seed = 84215

# DLR parameters
dlr_wmax = 10.0
dlr_eps = 1e-10

# === Run: Standard basis with anomalous chiAB ===
gf_struct = [('up', 1), ('dn', 1)]
S = Solver(beta=beta, gf_struct=gf_struct, dlr_wmax=dlr_wmax, dlr_eps=dlr_eps)

S.G0_iw['up'] << inverse(iOmega_n + mu)
S.G0_iw['dn'] << inverse(iOmega_n + mu)

h_int = U * n('up', 0) * n('dn', 0)

# Anomalous operators: Delta† = c†_up c†_dn, Delta = c_dn c_up
Delta_dag = c_dag('up', 0) * c_dag('dn', 0)
Delta = c('dn', 0) * c('up', 0)

S.solve(h_int=h_int,
        n_s=1,
        n_cycles=n_cycles,
        length_cycle=length_cycle,
        n_warmup_cycles=n_warmup,
        random_seed=seed,
        measure_chiAB_tau=True,
        chi_ops=[(Delta_dag, Delta)])

chi_tau = S.chiAB_tau

# === Exact result for single-site Hubbard atom ===
# States: |0> (E=0), |up> (E=-mu), |dn> (E=-mu), |ud> (E=-2mu+U)
# chi_pair(tau) = exp(E_d * (tau - beta)) / Z  for 0 < tau < beta
E0 = 0.0
E1 = -mu
Ed = -2 * mu + U
Z = np.exp(-beta * E0) + 2 * np.exp(-beta * E1) + np.exp(-beta * Ed)

# === Compare ===
if mpi.is_master_node():
    tau_mesh = chi_tau.mesh
    tau_vals = np.array([float(tau) for tau in tau_mesh])
    chi_data = chi_tau.data[:, 0].real

    # Exact values at the DLR tau points
    chi_exact = np.exp(Ed * (tau_vals - beta)) / Z

    print("=== Pairing susceptibility: anomalous chiAB vs exact ===")
    print(f"U = {U}, mu = {mu}, beta = {beta}")
    print(f"Z = {Z:.6f}, E_d = {Ed}")
    print(f"Number of DLR tau points: {len(tau_vals)}")
    print()
    print(f"{'tau':>10s} {'MC':>14s} {'exact':>14s} {'MC-exact':>14s}")
    for i, tau in enumerate(tau_vals):
        print(f"{tau:10.4f} {chi_data[i]:14.8f} {chi_exact[i]:14.8f} {chi_data[i] - chi_exact[i]:14.2e}")

    max_diff = np.max(np.abs(chi_data - chi_exact))
    print(f"\nMax |MC - exact|: {max_diff:.2e}")

    tol = 0.02
    assert max_diff < tol, f"Anomalous chiAB deviates from exact by {max_diff:.4f} > {tol}"
    print(f"PASSED: Anomalous chiAB agrees with exact result within tolerance {tol}.")
