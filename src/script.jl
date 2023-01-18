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
Rmax       = parse(Float64, ARGS[4]) # length associated to the solid
M          = parse(Int64,   ARGS[5]) # number of Chebyshev moments
NR         = parse(Int64,   ARGS[6]) # number of random vectors

# shapes = ["octahedron", "cube", "rhombic", "sphere"]
shape = "cube"

l = Rmax/(a_0/2) # number of atomic (100) planes

# Generate the positions of the atoms inside the nanoparticle
Elist, Edict, R = generate_shape(l, shape_name)

# Use those positions to generate the Hamiltonian matrix elements
H = gold_HV(Elist, Edict)

# Read the potential from file
Phi = comsol_read(potname)
