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
f.close()

Ni = len(Ei)
Nj = len(Ej)
print("Dimensions: ", Ni, Nj)

plt.pcolormesh(Ei,Ej,Φ)
plt.show()
