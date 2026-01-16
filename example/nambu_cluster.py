"""
CTINT reference calculation for 4x4 cluster with Nambu formalism (attractive U).

Nambu spinor: Psi = (c_up, c†_dn)^T
  - Psi_0 = c_up,   Psi†_0 = c†_up
  - Psi_1 = c†_dn,  Psi†_1 = c_dn

Nambu structure: single "nambu" block with 2 orbitals per site
  - orbital 2*i = particle sector for site i (Psi_0)
  - orbital 2*i+1 = hole sector for site i (Psi_1)

Physics: Anti-commutation relations require careful treatment
  n_up = c†_up c_up = Psi†_0 Psi_0 = n_0  (standard)
  n_dn = c†_dn c_dn = Psi_1 Psi†_1 = 1 - Psi†_1 Psi_1 = 1 - n_1

  U * n_up * n_dn = U * n_0 * (1 - n_1) = U * n_0 - U * n_0 * n_1

  Two terms:
  - Two-body: -U * n_0 * n_1 (sign flip from anti-commutation)
  - One-body: U * n_0 (absorbed in h_loc)

  h_loc structure:
  - particle sector: -mu + U
  - hole sector: +mu

  The solver interaction coefficient is -U (positive for attractive U < 0)

Correlators:
  G = -<T c_up(r, tau) c†_up(0, 0)>
  F = <T c_up(r, tau) c_dn(0, 0)> = -G_{particle, hole}
"""

from triqs_ctint import Solver
from triqs.gf import *
from triqs.gf.meshes import MeshCycLat
from triqs.lattice import BravaisLattice
from triqs.operators import n
from h5 import HDFArchive
import numpy as np
import triqs.utility.mpi as mpi

# Physical parameters - attractive Hubbard model
U = -2.0   # Attractive interaction
t = 1.0    # Nearest-neighbor hopping (H = -t Σ c†c)
mu = U / 2 # Chemical potential (half-filling)
beta = 5.0

# DLR parameters
w_max = 50.0
eps = 1e-10

# Lattice geometry
Nx, Ny = 4, 4
n_sites = Nx * Ny

# Hopping matrix with periodic boundary conditions (n_sites x n_sites)
hopping = np.zeros((Nx, Ny, Nx, Ny))
for d in range(Nx):
    hopping[d, :, (d + 1) % Nx, :] = t * np.eye(Ny)
    hopping[(d + 1) % Nx, :, d, :] = t * np.eye(Ny)
for d in range(Ny):
    hopping[:, d, :, (d + 1) % Ny] = t * np.eye(Nx)
    hopping[:, (d + 1) % Ny, :, d] = t * np.eye(Nx)
hopping = hopping.reshape(n_sites, n_sites)

# Nambu structure: single block with 2 orbitals (particle=0, hole=1) per site
n_orb = 2 * n_sites

# Build Nambu h_loc matrix (n_orb x n_orb)
# hopping[i,j] = t for neighbors, so -hopping gives the kinetic term H = -t Σ c†c
# h_loc structure per site:
#   particle (2*i):   -mu + U - hopping[i,j]
#   hole (2*i+1):     +mu + hopping[i,j]  (flipped dispersion)
hloc = np.zeros((n_orb, n_orb))
for i in range(n_sites):
    for j in range(n_sites):
        # particle-particle block (H = -t Σ c†c)
        hloc[2*i, 2*j] = -hopping[i, j]
        if i == j:
            hloc[2*i, 2*i] += -mu + U
        # hole-hole block (flipped sign for dispersion)
        hloc[2*i+1, 2*j+1] = +hopping[i, j]
        if i == j:
            hloc[2*i+1, 2*i+1] += +mu

# Interaction Hamiltonian in Nambu basis: -U * n_0 * n_1
# Solver expands in powers of (-U_input), so U_input = -U
h_int = sum(-U * n("nambu", 2*i) * n("nambu", 2*i+1) for i in range(n_sites))

# Block structure: single nambu block
gf_struct = [("nambu", n_orb)]

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=100, n_tau=201)

# Initialize non-interacting Green's function
S.G0_iw["nambu"] << inverse(iOmega_n - hloc)

# Solve
S.solve(h_int=h_int,
        delta=0.1,
        n_s=1, 
        n_cycles=100000,
        length_cycle=100,
        n_warmup_cycles=1000,
        measure_histogram=True,
        measure_M_tau=False,
        measure_M_iw_dlr=True,
        w_max=w_max,
        eps=eps,
        measure_density=True)

# Convert G_iw_dlr to DLR imaginary time
G_dlr_tau = make_gf_dlr_imtime(S.G_iw_dlr)

# Extract tau mesh
tau_mesh = G_dlr_tau.mesh

# Lattice mesh for output
clat_mesh = MeshCycLat(BravaisLattice([[1, 0], [0, 1]]), [Nx, Ny, 1])

# Helper: site index to (ix, iy)
def site_to_xy(i):
    return i // Ny, i % Ny

# Helper: (rx, ry) to linear r_idx in clat_mesh ordering
def xy_to_ridx(rx, ry):
    return rx * Ny + ry

# Translational averaging for normal Green's function
# G(r, tau) = (1/N) sum_i G_{particle(i), particle(i+r)}(tau)
def extract_G_r_tau(G_block_tau):
    """Extract G(r, tau) from particle-particle block with translational averaging."""
    result = Gf(mesh=MeshProduct(clat_mesh, tau_mesh), target_shape=[])

    for rx in range(Nx):
        for ry in range(Ny):
            r_idx = xy_to_ridx(rx, ry)
            # Average over all sites i
            avg = np.zeros(len(tau_mesh), dtype=complex)
            for i in range(n_sites):
                ix, iy = site_to_xy(i)
                jx, jy = (ix + rx) % Nx, (iy + ry) % Ny
                j = jx * Ny + jy
                # particle-particle: orbital 2*i, 2*j
                avg += G_block_tau.data[:, 2*i, 2*j]
            result.data[r_idx, :] = avg / n_sites

    return result

# Translational averaging for anomalous correlator
# F(r, tau) = <T c_up(r,tau) c_dn(0,0)> = -G_{particle, hole}
# (TRIQS convention: G_ab = -<T Psi_a Psi†_b>, so G_01 = -<T c_up c_dn> = -F)
def extract_F_r_tau(G_block_tau):
    """Extract F(r, tau) from particle-hole block with translational averaging."""
    result = Gf(mesh=MeshProduct(clat_mesh, tau_mesh), target_shape=[])

    for rx in range(Nx):
        for ry in range(Ny):
            r_idx = xy_to_ridx(rx, ry)
            # Average over all sites i
            avg = np.zeros(len(tau_mesh), dtype=complex)
            for i in range(n_sites):
                ix, iy = site_to_xy(i)
                jx, jy = (ix + rx) % Nx, (iy + ry) % Ny
                j = jx * Ny + jy
                # particle-hole: orbital 2*i, 2*j+1
                # F = -G_01 (sign from TRIQS convention)
                avg -= G_block_tau.data[:, 2*i, 2*j+1]
            result.data[r_idx, :] = avg / n_sites

    return result

# Extract G and F from the single nambu block
G_nambu_tau = G_dlr_tau["nambu"]
G_r_tau = extract_G_r_tau(G_nambu_tau)
F_r_tau = extract_F_r_tau(G_nambu_tau)

n_r = len(clat_mesh)
n_tau = len(tau_mesh)

# Print some results
if mpi.is_master_node():

    print(f"\n=== Nambu Green's Functions (Attractive Hubbard) ===")
    print(f"U = {U}, t = {t}, mu = {mu}, beta = {beta}")
    print(f"Lattice: {Nx}x{Ny}, n_r = {n_r}, n_tau = {n_tau}")

    # Find tau closest to beta/2
    tau_vals = [float(tau_mesh.to_value(i)) for i in range(n_tau)]
    mid_idx = min(range(len(tau_vals)), key=lambda i: abs(tau_vals[i] - beta/2))
    tau_mid = tau_vals[mid_idx]

    print(f"\nNormal Green function G(R=0, tau) samples:")
    for tau_idx in range(min(n_tau, 10)):
        tau = tau_vals[tau_idx]
        print(f"  tau={tau:8.4f}: G={G_r_tau.data[0, tau_idx].real:12.6f}")

    print(f"\nAnomalous correlator F(R=0, tau) samples:")
    for tau_idx in range(min(n_tau, 10)):
        tau = tau_vals[tau_idx]
        print(f"  tau={tau:8.4f}: F={F_r_tau.data[0, tau_idx].real:12.6f}")

    print(f"\nG(R, tau={tau_mid:.4f}) samples:")
    for r_idx in range(min(n_r, 10)):
        rx, ry = r_idx // Ny, r_idx % Ny
        print(f"  R=({rx},{ry},0): G={G_r_tau.data[r_idx, mid_idx].real:12.6f}")

    print(f"\nF(R, tau={tau_mid:.4f}) samples:")
    for r_idx in range(min(n_r, 10)):
        rx, ry = r_idx // Ny, r_idx % Ny
        print(f"  R=({rx},{ry},0): F={F_r_tau.data[r_idx, mid_idx].real:12.6f}")

# Save results
if mpi.is_master_node():
    print(f"\nWriting results to nambu_cluster.h5 ...")
    with HDFArchive("nambu_cluster.h5", 'w') as arch:
        arch["U"] = U
        arch["t"] = t
        arch["mu"] = mu
        arch["beta"] = beta
        arch["n_r"] = n_r
        arch["n_tau"] = n_tau
        arch["S"] = S
        arch["G_r_tau"] = G_r_tau
        arch["F_r_tau"] = F_r_tau
