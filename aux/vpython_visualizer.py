import vpython as vp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import sys

# Load the file with the atomic positions. If it has a potential, 
# it has to be treated with complex numbers and has 4 columns
hasV = False
filename = sys.argv[1]

# Fetch the data from the comsol file
positions = []
Vset = []
with open(filename, 'r') as f:
    a = f.readlines() 

    b = a[9:] # remove the header

    # check the size of the first data line. If it has 4 entries, it has a potential
    if len(b[0].split()) == 4:
        hasV = True

    for line in b:
        listed = line.split()
        x = float(listed[0])
        y = float(listed[1])
        z = float(listed[2])
        V = listed[-1]
        positions.append([x,y,z])
        Vset.append(V)

N = len(positions)
print("Number of atoms: ", N)
if hasV: print("The complex potential is specified")
else:    print("The complex potential is not specified")

# Convert the data into a numpy array
positions = np.array(positions)

# If the potential exists, convert it into an array
Varray = np.zeros(N, dtype=complex)
if hasV:
    for i,v in enumerate(Vset):
        c = v.split(".")
        re = c[0] + "." + c[1][:-2]
        im = c[1][-2:] + "." +  c[2][:-1]
        V = float(re) + 1j*float(im)
        Varray[i] = V

# Represent the average complex potential along the z axis
if hasV:

    # find all unique zs and index them with a dictionary
    zset = list(set(positions[:,2]))
    zdic = {z:i for i,z in enumerate(zset)}
    Nzs  = len(zset)

    # iterate over every atomic position and project it onto z
    vz = np.zeros(Nzs, dtype=complex)
    Nsites = np.zeros(Nzs, dtype=int)
    for i in range(N):
        z = positions[i,2]
        V = Varray[i]
        j = zdic[z]
        vz[j] += V
        Nsites[j] += 1

    vz = vz/Nsites # perform the average

    # order the sets
    inds = sorted(range(Nzs), key=lambda k: float(zset[k]))
    zset = np.array(zset)[inds]
    vz = vz[inds]

    # plot
    fig, axs = plt.subplots(1,1,figsize=(8,5))
    axs.plot(zset, np.real(vz), label="Real")
    axs.plot(zset, np.imag(vz), label="Imag")
    axs.set_title("Average complex potential along z")
    axs.legend()
    plt.show()




# Select whether you want the real or the imaginary part
Vs = np.real(Varray)
# Vs = np.imag(Varray)

# Use the complex potential to create the colors. If no potential exists, the color is white
colors = [[1,1,1] for i in range(N)]
if hasV:
    maxV = np.max(Vs)
    minV = np.min(Vs)
    cmap = cm.get_cmap('viridis',12)
    colors = cmap((Vs-minV)/(maxV-minV))[:,:3]


vp.canvas(title='Nanoparticle', width=1200, height=800)

a = 2.0*np.max(positions)
radius = 0.6*a/N**0.33
w = radius/10.0
ex = vp.arrow(axis=vp.vector(a,0,0),color=vp.color.red, shaftwidth=w)
ey = vp.arrow(axis=vp.vector(0,a,0),color=vp.color.green, shaftwidth=w)
ez = vp.arrow(axis=vp.vector(0,0,a),color=vp.color.blue, shaftwidth=w)

# display all the atoms with the color given by the potential
for col,pos in zip(colors,positions):

    vec = vp.vector(*pos)
    col = vp.vec(*col)

    vp.sphere(pos=vec, radius=radius, color=col)

