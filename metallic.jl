# redundancy defining xx and yy in aux.jl while summing the moments

using DataStructures
using Interpolations
using HDF5
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
HaeV = 27.211386245988       # [eV]        Hartree in eV
boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant
#a_0 = 7.291                  # Lattice constant

# Important parameters
Fermi = 0.2595               # [Ha]     Fermi energy
T     = 298                  # [Kelvin] room temperature 
beta  = 1/(boltzmann*T/HaeV) # [Ha⁻¹]   inverse temperature

# Modification to the spectral lines
kappa = 0.0
gph   = 0.06/HaeV            # [Ha] broadening

# List of frequencies

hwlist=LinRange(2.0,4.0,21)|>collect
hwlist=hwlist./HaeV
PBS_INDEX=get(ENV,"PBS_A","1")
PBS_INDEX=parse(Int64,PBS_INDEX)

# these are exactly the same as hwlist, but in string format
filelist=["2p0eV","2p1eV","2p2eV","2p3eV","2p4eV","2p5eV","2p6eV",
          "2p7eV","2p8eV","2p9eV","3p0eV","3p1eV","3p2eV","3p3eV",
          "3p4eV","3p5eV","3p6eV","3p7eV","3p8eV","3p9eV","4p0eV",
          "4p1eV","4p2eV","4p3eV","4p4eV","4p5eV","4p6eV","4p7eV",
          "4p8eV","4p9eV"]
hw=hwlist[PBS_INDEX]
file_hw = filelist[PBS_INDEX]

include("aux.jl")
include("shape_lib.jl")
include("potential.jl")
include("gold.jl")

# parse input
Nargs = size(ARGS)[1]
println("Nargs:",Nargs)
shape_name = ARGS[1]
Rmax       = parse(Float64, ARGS[2])
hw         = parse(Float64, ARGS[3])/HaeV
outfile_name = ARGS[4]
h5name = ARGS[5]
println("hw",hw)

# Rmax = 2.0 # 'radius' of the nanoparticle. If it's not a sphere, it controls the size
# shape_name = "cube"
# shape_name = "octahedron"
# shape_name = "rhombic_dodecahedron"

# Generate the positions of the atoms inside the nanoparticle
l = Rmax/(a_0/2) # number of atomic (100) planes
Elist, Edict, R = generate_shape(l, shape_name)
println("Finished creating the atomic positions. Number of atoms: ", length(Elist))
# print the positions of the particles for visualization
# print_positions(R)

# Use those positions to generate the Hamiltonian matrix elements
H = gold_HV(Elist, Edict)

Phi = compute_potential_sphere(R, Edict)
# If the potential file is provided, use it instead of the default linear potential
if Nargs == 6
    potname = ARGS[6]
    Phi = comsol_read(potname)
    # Make sure that this potential is compatible with the system 
    if size(Phi)!=size(H)
        println("Dimensions of the potential provided do not match the geometry requested. Exiting.")
        exit()
    end
end


# a,b=size(H)
# println("Dimensions of the sparse Hamiltonian: $a $b")
#exact(H, Phi, Fermi, beta)


M1 = 30 # number of Cheb moments PER BLOCK for first operator 
M2 = 30 # number of Cheb moments PER BLOCK for second operator
NR = 1 # number of random vectors
D1 = 40 # number of row blocks in the mu matrix
D2 = 40 # number of col blocks in the mu matrix

println("Computing Chebyshev moments")
mumn=compute_mumn!(H,Phi,M1,M2,D1,D2,NR)
# writedlm("cheb_moments.dat", mumn)

# h5name = "moments.h5"
fid = h5open(h5name, "w")
fid["moments"] = mumn
close(fid)

# h5name = "moments.h5"
# fid = h5open(h5name, "r")
# mumn = read(fid["moments"])
# close(fid)

# N1 = 1000
# N2 = N1

# println("Resumming Chebyshev moments")
# Ej,Nelist,Nhlist=compute_eh_list(mumn,N1,N2,hw,goldA,goldB,Fermi,beta,kappa,gph)

# f3=open(outfile_name,"w")
# for i=1:lastindex(Ej)
    # write(f3,@sprintf(" %+1.15E",Ej[i]))
    # write(f3,@sprintf(" %+1.15E",Nelist[i]))
    # write(f3,@sprintf(" %+1.15E",Nhlist[i]))
    # write(f3,"\n")
# end


