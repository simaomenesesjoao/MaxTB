# redundancy defining xx and yy in aux.jl while summing the moments

using DataStructures
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
# boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant
# a_0 = 0.408 # Lattice constant (length of FCC unit cell in nanometers)

include("aux.jl")
include("shape_lib.jl")
include("potential.jl")
include("gold.jl")

# parse input
Nargs = size(ARGS)[1]
if Nargs != 6
    println("Number of arguments (", Nargs, ") is wrong. Usage:\n")
    println("julia moments.jl [potential filename] [output h5df filename] [shape] [length] [M] [NR].\n")
    println("Example: julia moments.jl potential.dat cheb_moments.h5 sphere 2.0 1000 3\n")
    println("will calculate the Chebyshev matrix with 1000 polynomials and 3 random vectors for a sphere of radius 2nm. The potential will be read from file 'potential.dat' and the Chebyshev matrix will be written to file 'cheb_moments.h5'.")
    exit()
end

potname    = ARGS[1] # name of the input potential file
h5name     = ARGS[2] # name of the output hdf5 file with the Chebyshev moments
shape_name = ARGS[3] # name of the shape
Rmax       = parse(Float64, ARGS[4]) # length associated to the solid
M          = parse(Int64,   ARGS[5]) # number of Chebyshev moments
NR         = parse(Int64,   ARGS[6]) # number of random vectors

shapes = ["octahedron", "cube", "rhombic", "sphere"]
if !(shape_name in shapes)
    println("Shape "*shape_name*" not supported. Use sphere, octahedron, cube or rhombic")
    exit()
end

# Generate the positions of the atoms inside the nanoparticle
l = Rmax/(a_0/2) # number of atomic (100) planes
Elist, Edict, R = generate_shape(l, shape_name)
# R = zeros(Float64,NN,3)
# N,b = size(R)
# println("shape:",N," ",b)
# println("Finished creating the atomic positions. Number of atoms: ", length(Elist))

# println("Edict:\n", Edict)
# println("Elist:\n", Elist)
# println("R:\n", R)

# Edict2 = Dict{Vector{Int64},Int64}()
# Elist2 = Deque{Vector{Int64}}()
# pos_int = [0, 0, 0]
# for i in 1:N
    # pos = R[i,:]
    # pos_int[1] = round(pos[1]*2/a_0)
    # pos_int[2] = round(pos[2]*2/a_0)
    # pos_int[3] = round(pos[3]*2/a_0)
    
    # Edict2[pos_int] = i
    # push!(Elist2, pos_int)
    # println("position: ", i, " ", pos, pos_int)
    # println("new Edict:", Edict2[pos_int])
    # println("original Edict: ", Edict[pos_int])
    # println("\n")
# end

# println("Compare Elists")
# println("New Elist:\n",Elist2,"\n")
# println("Old Elist:\n",Elist,"\n")


# Use those positions to generate the Hamiltonian matrix elements
H = gold_HV(Elist, Edict)

# Read the potential from file
Phi = comsol_read(potname)

# Make sure that this potential is compatible with the system 
a,b=size(H)
println("Dimensions of the sparse Hamiltonian: $a $b")
if size(Phi)!=size(H)
    println("Dimensions of the potential provided do not match the geometry requested. Exiting.")
    exit()
end

#exact(H, Phi, Fermi, beta)

M1 = 30 # number of Cheb moments PER BLOCK for first operator 
M2 = 30 # number of Cheb moments PER BLOCK for second operator
D1 = ceil(Int, M/M1)
D2 = ceil(Int, M/M2)
Meff = D1*M1

println("Computing Chebyshev moments")
println("Computing ",Meff,"x",Meff," Chebyshev matrix in blocks of size ",M1,"x",M2," for a total of ",D1,"x",D2," blocks.")
mumn = compute_mumn!(H,Phi,M1,M2,D1,D2,NR)

println("Storing Chebyshev moments in hdf file: ", h5name)
fid = h5open(h5name, "w")
fid["moments"] = mumn
close(fid)
