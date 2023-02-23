import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import sys

# file = "/mnt/data/projects_sync/imperial_nanoparticles/nanoparticle-fgr/examples/phi.h5"
file = sys.argv[1]

f = h5.File(file, 'r')
Ei = f["Ei"][:]
Ej = f["Ej"][:]
Φ  = f["Phi"][:,:]

xs = f["xs"][:]
ys = f["ys"][:]

Ke = f["Ke"][:]
Kh = f["Kh"][:]
Phi_itp = f["Phi_itp"][:]
mask_el = f["mask_el"][:]
mask_ho = f["mask_ho"][:]

Ne = f["Nelist"][:]
Nh = f["Nhlist"][:]

f.close()

# Ni = len(Ei)
# Nj = len(Ej)
print("plottning: ")

fig, axs = plt.subplots(2,3,figsize=(10,5))
axs[0,0].pcolormesh(Ei, Ej, Φ, rasterized=True)
axs[1,0].pcolormesh(xs, ys, Phi_itp, rasterized=True)
axs[0,1].pcolormesh(xs, ys, mask_el, rasterized=True)
axs[1,1].pcolormesh(xs, ys, Ke, rasterized=True)
axs[0,2].pcolormesh(xs, ys, mask_ho, rasterized=True)
axs[1,2].pcolormesh(xs, ys, Kh, rasterized=True)

plt.show()
