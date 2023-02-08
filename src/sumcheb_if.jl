# Interface to the sumcheb.jl code
using HDF5
using DelimitedFiles
using Printf

include("sumcheb.jl")


Nargs = size(ARGS)[1]
if Nargs != 8
    println("Number of arguments (", Nargs, ") is wrong. Usage:\n")
    println("julia sumcheb_if.jl [input hdf filename] [output hdf filename] [A] [B] [NE1] [NE2] [percentage 1] [percentage 2].\n")
    println("Example: julia process_moments.jl hdfile.h5 phi.dat 0.2 0.7 1000 1100 70 80\n")
    println("will calculate optical matrix in energy space (phi.dat) using the Chebyshev moments in the hdf file (hdffile.h5) with 1000 energy points in the first index and 1100 energy points in the second index. 70% of the polynomials in the first index will be used and 80% of the polynomials in the second index will be used. The Hamiltonian's scaling parameters are A (shift) and B (scale). The output is one single file with the list of energies in the first index, the list of energies in the second index and then the optical matrix.")
    exit()
end

h5name  = ARGS[1] # name of the file with the Chebyshev moments
outname = ARGS[2] # name of the output file with the optical matrix
A       = parse(Float64, ARGS[3]) # Hamiltonian scaling factor A (shift)
B       = parse(Float64, ARGS[4]) # Hamiltonian scaling factor B (scale)
N1      = parse(Int,     ARGS[5]) # Number of energy points in first index
N2      = parse(Int,     ARGS[6]) # Number of energy points in second index
perc1   = parse(Float64, ARGS[7]) # percentage of the Cheb polys in first index
perc2   = parse(Float64, ARGS[8]) # percentage of the Cheb polys in second index

println("Arguments read from command line: ", h5name, " ", outname, " ", A, " ", B, " ", N1, " ", N2, " ", perc1, " ", perc2, "\n")


# read from the hdf5 file
fid = h5open(h5name, "r")
mumn = read(fid["moments"])
close(fid)

println("Resumming Chebyshev moments")
Ei, Ej, Phi = resum_mu_dd(mumn, N1, N2, A, B, perc1, perc2)

println("Storing in hdf file: ", h5name)
fid = h5open(outname, "w")
fid["Phi"] = Phi
# fid["Phi"] = 1 .+ 0 .* Phi
fid["Ei"]  = Ei |> collect
fid["Ej"]  = Ej |> collect
close(fid)

# Testing norm
dif = test_conv(Phi, Ei, Ej)
println("norm: ", dif)
