import vpython as vp
import numpy as np
from matplotlib import cm
import sys

# Load the file with the atomic positions. If it has a potential, 
# it has to be trated with complex numbers and has 4 columns
hasV = True
filename = "/mnt/data/projects_sync/imperial_nanoparticles/" + sys.argv[1]
atoms = np.loadtxt(filename, dtype=complex)
if len(atoms[0]) == 3:
    hasV = False
    atoms = np.loadtxt(filename)

N = len(atoms)
print("Number of atoms: ", N)

xs = np.real(atoms[:,0])
ys = np.real(atoms[:,1])
zs = np.real(atoms[:,2])
Vs = np.real(atoms[:,-1])

maxV = np.max(Vs)
minV = np.min(Vs)

cmap = cm.get_cmap('viridis',12)
colors = cmap((Vs-minV)/(maxV-minV))[:,:3]

vp.canvas(title='Nanoparticle', width=1200, height=800)

a = 2.0*np.max(xs)
radius = 0.6*a/N**0.33
w = radius/10.0
ex = vp.arrow(axis=vp.vector(a,0,0),color=vp.color.red, shaftwidth=w)
ey = vp.arrow(axis=vp.vector(0,a,0),color=vp.color.green, shaftwidth=w)
ez = vp.arrow(axis=vp.vector(0,0,a),color=vp.color.blue, shaftwidth=w)

# display all the atoms
for i,pos in enumerate(atoms):
    x = xs[i]
    y = ys[i]
    z = zs[i]
    col = colors[i]

    vec = vp.vector(x,y,z)
    col = vp.vec(*col)
    if not hasV:
        col = vp.color.white


    vp.sphere(pos=vec, radius=radius, color=col)

