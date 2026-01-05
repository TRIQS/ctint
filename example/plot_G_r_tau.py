from triqs.gf import *
from h5 import HDFArchive
import matplotlib.pyplot as plt
import numpy as np

with HDFArchive("cluster.h5", 'r') as arch:
    U = arch["U"]
    G_r_tau = arch["G_r_tau"]
    G0_r_tau = arch["G0_r_tau"]

# Extract G(r=0, tau) on DLR mesh
dlr_tau_mesh = G_r_tau.mesh.components[1]
dlr_tau_vals = np.array([t.value for t in dlr_tau_mesh])
G_local_dlr = G_r_tau.data[0, :]
G0_local_dlr = G0_r_tau.data[0, :]

# Interpolate to dense uniform tau grid
def to_dense_imtime(data, mesh, n_tau=1001):
    gf = Gf(mesh=mesh, target_shape=[])
    gf.data[:] = data
    return make_gf_imtime(make_gf_dlr(gf), n_tau)

G_local_dense = to_dense_imtime(G_local_dlr, dlr_tau_mesh)
G0_local_dense = to_dense_imtime(G0_local_dlr, dlr_tau_mesh)
dense_tau_vals = np.array([t.value for t in G_local_dense.mesh])

plt.figure(figsize=(8, 5))
plt.plot(dense_tau_vals, G_local_dense.data.real, '-', label=r'$G(r=0, \tau)$')
plt.plot(dlr_tau_vals, G_local_dlr.real, 'o', markersize=6, label='DLR grid points')
plt.plot(dense_tau_vals, G0_local_dense.data.real, '--', label=r'$G_0(r=0, \tau)$ ($\mu=0$)')
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G(\tau)$')
plt.title(f'Local Green\'s function ($U = {U}$)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('G_r0_tau.pdf')
