using HDF5
using DelimitedFiles
using Printf

include("dist.jl")
HaeV = 27.211386245988       # [eV]        Hartree in eV

Nargs = size(ARGS)[1]
if Nargs != 8
    println("Number of arguments is 8")
    println(ARGS)
end
h5name  = ARGS[1] # name of the file with the optical matrix
outname = ARGS[2] # name of the output file with the electron distribution
N1      = parse(Int,     ARGS[3]) # Number of energy points in first index
N2      = parse(Int,     ARGS[4]) # Number of energy points in second index
fermi   = parse(Float64, ARGS[5]) # Fermi energy
hw      = parse(Float64, ARGS[6]) # Frequency
kbT     = parse(Float64, ARGS[7]) # Temperature
debug   = parse(Int,     ARGS[8]) # Debug flag to print out the energy matrices

# read from the hdf5 file
fid = h5open(h5name, "r")
Phi = read(fid["Phi"])
Ei  = read(fid["Ei"])*HaeV
Ej  = read(fid["Ej"])*HaeV
close(fid)

Gammaf = 1

# Test different temperatures and numbers of integration points
Nrows = 1
if debug > 0
    Nrows = 5
end

columns = Array{Float64}(undef, N1, Nrows)
columns[:,1] = LinRange(Ej[1], Ej[end], N1)
columns[:,2] = sum_dist(Phi, Ei, Ej, N1, N2,   fermi, hw, kbT,   Gammaf, debug)

if debug > 0
    columns[:,3] = sum_dist(Phi, Ei, Ej, N1, N2*3, fermi, hw, kbT,   Gammaf, false)
    columns[:,4] = sum_dist(Phi, Ei, Ej, N1, N2,   fermi, hw, kbT*2, Gammaf, false)
    columns[:,5] = sum_dist(Phi, Ei, Ej, N1, N2,   fermi, hw, kbT,   Gammaf*2, false)
end

# Write the data
writedlm(outname, columns)
