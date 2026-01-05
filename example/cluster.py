from triqs_ctint import Solver
from triqs.gf import *
from triqs.gf.meshes import MeshCycLat
from triqs.lattice import BravaisLattice
from triqs.operators import n
from h5 import HDFArchive
import numpy as np
import triqs.utility.mpi as mpi

# Physical parameters
U = 2.0  # Density-density interaction
t = 1.0  # Hopping
mu = U/2.0 + 0.4  # Chemical Potential
beta = 5.0  # Inverse temperature

# DLR parameters
w_max = 50.0
eps = 1e-10

# Lattice geometry
Nx, Ny = 4, 4
n_orb = Nx * Ny

# Hopping matrix with periodic boundary conditions
hopping = np.zeros((Nx, Ny, Nx, Ny))
for d in range(Nx):
    hopping[d, :, (d + 1) % Nx, :] = t * np.eye(Ny)
    hopping[(d + 1) % Nx, :, d, :] = t * np.eye(Ny)
for d in range(Ny):
    hopping[:, d, :, (d + 1) % Ny] = t * np.eye(Nx)
    hopping[:, (d + 1) % Ny, :, d] = t * np.eye(Nx)
hopping = hopping.reshape(n_orb, n_orb)

hloc0 = -mu * np.eye(n_orb) - hopping
hloc0_mu0 = -hopping  # For G0 reference (mu=0)

# Interaction Hamiltonian
h_int = sum(U * n("up", i) * n("dn", i) for i in range(n_orb))

# Density operators for chi measurement (up spin only)
n_ops = [n("up", i) for i in range(n_orb)]

# Block structure
gf_struct = [("dn", n_orb), ("up", n_orb)]

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=100, n_tau=201)

# Initialize non-interacting Green's function
for bl, g_bl in S.G0_iw:
    g_bl << inverse(iOmega_n - hloc0)

# Solve using M_iw_dlr measurement
S.solve(h_int=h_int,
        delta=0.1,
        n_cycles=100000,
        length_cycle=100,
        n_warmup_cycles=1000,
        measure_histogram=True,
        measure_M_tau=False,
        measure_M_iw_dlr=True,
        w_max=w_max,
        eps=eps,
        measure_density=True,
        measure_chiAB_tau=True,
        chi_A_vec=n_ops,
        chi_B_vec=n_ops,
        n_tau_chi2=201)

# Convert G_iw_dlr to DLR imaginary time
G_dlr_tau = make_gf_dlr_imtime(S.G_iw_dlr)

# Compute G0(r, tau) with mu=0 for reference
G0_iw_dlr = BlockGf(mesh=S.G_iw_dlr.mesh, gf_struct=gf_struct)
for bl, g_bl in G0_iw_dlr:
    g_bl << inverse(iOmega_n - hloc0_mu0)
G0_dlr_tau = make_gf_dlr_imtime(G0_iw_dlr)

# Translational averaging: F(r) = (1/N) * sum_i F_{i, i+r}
clat_mesh = MeshCycLat(BravaisLattice([[1, 0], [0, 1]]), [Nx, Ny, 1])
ix, iy = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')

def translational_average(data_arr, tau_mesh):
    """Average F_{ij}(tau) over translations: F(r, tau) = (1/N) sum_i F_{i,i+r}(tau)"""
    result = Gf(mesh=MeshProduct(clat_mesh, tau_mesh), target_shape=[])
    data_4d = data_arr.reshape(len(tau_mesh), Nx, Ny, Nx, Ny)
    for rx in range(Nx):
        for ry in range(Ny):
            jx, jy = (ix + rx) % Nx, (iy + ry) % Ny
            result.data[rx * Ny + ry, :] = data_4d[:, ix, iy, jx, jy].sum(axis=(1, 2)) / n_orb
    return result

# G(r, tau) averaged over blocks
tau_mesh = G_dlr_tau.mesh
G_r_tau = translational_average(sum(g.data for _, g in G_dlr_tau), tau_mesh)
G_r_tau.data[:] /= len(gf_struct)

G0_r_tau = translational_average(sum(g.data for _, g in G0_dlr_tau), tau_mesh)
G0_r_tau.data[:] /= len(gf_struct)

# Connected chi(r, tau) = <n_up(r,tau) n_up(0,0)> - <n_up>^2
chi_r_tau = translational_average(S.chiAB_tau.data, S.chiAB_tau.mesh)
n_up_avg = S.density[1].trace() / n_orb  # index 1 = "up" block
chi_r_tau.data[:] -= n_up_avg**2

# Save results
if mpi.is_master_node():
    with HDFArchive("cluster.h5", 'w') as arch:
        arch["U"] = U
        arch["S"] = S
        arch["G_r_tau"] = G_r_tau
        arch["G0_r_tau"] = G0_r_tau
        arch["chi_r_tau"] = chi_r_tau
