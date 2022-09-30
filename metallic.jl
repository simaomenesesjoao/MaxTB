# redundancy defining xx and yy in aux.jl while summing the moments

using DataStructures
using Interpolations
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using MKLSparse
using Dates
using Random
using MKL
using Printf
using NearestNeighbors


# Important constants
HaeV = 27.211386245988    # Hartree in eV
boltzmann = 8.617333e-5   # Boltzmann constant in eV/Kelvin
a_0 = 7.291               # Lattice constant

# Important parameters
Fermi=0.5190/2            # Fermi energy
T = 298                   # room temperature 
beta=1/(boltzmann*T/HaeV) # inverse temperature in Hartree


# List of frequencies
hwlist=LinRange(2.0,4.0,21)|>collect
hwlist=hwlist./HaeV
PBS_INDEX=get(ENV,"PBS_ARRAY_INDEX","1")
PBS_INDEX=parse(Int64,PBS_INDEX)

# these are exactly the same as hwlist, but in string format
filelist=["2p0eV.txt","2p1eV.txt","2p2eV.txt","2p3eV.txt","2p4eV.txt","2p5eV.txt","2p6eV.txt",
          "2p7eV.txt","2p8eV.txt","2p9eV.txt","3p0eV.txt","3p1eV.txt","3p2eV.txt","3p3eV.txt",
          "3p4eV.txt","3p5eV.txt","3p6eV.txt","3p7eV.txt","3p8eV.txt","3p9eV.txt","4p0eV.txt",
          "4p1eV.txt","4p2eV.txt","4p3eV.txt","4p4eV.txt","4p5eV.txt","4p6eV.txt","4p7eV.txt",
          "4p8eV.txt","4p9eV.txt"]
hw=hwlist[PBS_INDEX]
file_hw = filelist[PBS_INDEX]

include("aux.jl")
include("generate_positions.jl")
include("potential.jl")
include("gold.jl")

Rmax = 4.0 # 'radius' of the nanoparticle. If it's not a sphere, it controls the size
# shape_name = "cube"
shape_name = "octahedron"
# shape_name = "rhombic_dodecahedron"

# Generate the positions of the atoms inside the nanoparticle
Elist, Edict, R = generate_shape(Rmax, shape_name)
println("Finished creating the atomic positions. Number of atoms: ", length(Elist))
# print the positions of the particles for visualization
# print_positions(R)

# Use those positions to generate the Hamiltonian matrix elements
H = gold_HV(Elist, Edict)

Phi = compute_potential_sphere(R, Edict)
# Phi,PhiT = comsol_read("comsol_out.txt")


# a,b=size(H)
# println("Dimensions of the sparse Hamiltonian: $a $b")
# exact(H, Phi, Fermi, beta)


M1 = 200 # number of Cheb moments for first operator
M2 = 200 # number of Cheb moments for second operator
NR = 5 # number of random vectors
D1 = 1 # number of row blocks in the mu matrix
D2 = 1 # number of col blocks in the mu matrix

println("Computing Chebyshev moments")
mumn=compute_mumn!(H,Phi,M1,M2,D1,D2,NR)
writedlm("cheb_moments.dat", mumn)

kappa = 0.0
gph = 0.06/HaeV
N1 = 1000
N2 = N1

println("Resumming Chebyshev moments")
Ej,Nelist,Nhlist=compute_eh_list(mumn,N1,N2,hw,goldA,goldB,Fermi,beta,kappa,gph)

f3=open("Metallic_output_"*file_hw*".txt","w")
for i=1:lastindex(Ej)
    write(f3,@sprintf(" %+1.15E",Ej[i]))
    write(f3,@sprintf(" %+1.15E",Nelist[i]))
    write(f3,@sprintf(" %+1.15E",Nhlist[i]))
    write(f3,"\n")
end


