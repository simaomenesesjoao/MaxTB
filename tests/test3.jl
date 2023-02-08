# Testing the code, using julia scripts instead of bash

using HDF5
using DelimitedFiles

include("../src/sumcheb.jl")
include("../src/dist.jl")
include("../src/hcg.jl")

HaeV = 27.211386245988

# Parameters to run the first part of the script: converting moments into 
# optical matrix in energy space
h5name = "moments.h5"
outname = "phi.h5"
A = 0.1
B = 1.0
N1 = 1024
N2 = 1024
perc1 = 100
perc2 = 70
perc3 = 50

fid = h5open(h5name, "r")
mumn = read(fid["moments"])
close(fid)

println("Resumming Chebyshev moments")
Ei, Ej, Phi  = resum_mu_dd(mumn, N1, N2, A, B, perc1, perc1)

println("Resumming Chebyshev moments perc2")
Ei, Ej, Phi2 = resum_mu_dd(mumn, N1, N2, A, B, perc2, perc2)

println("Resumming Chebyshev moments perc3")
Ei, Ej, Phi3 = resum_mu_dd(mumn, N1, N2, A, B, perc3, perc3)

Ei *= HaeV
Ej *= HaeV

fid = h5open("Phi.h5", "w")
fid["Phi"] = Phi
fid["Ei"]  = Ei |> collect
fid["Ej"]  = Ej |> collect
close(fid)

fid = h5open("Phi2.h5", "w")
fid["Phi"] = Phi2
fid["Ei"]  = Ei |> collect
fid["Ej"]  = Ej |> collect
close(fid)

fid = h5open("Phi3.h5", "w")
fid["Phi"] = Phi3
fid["Ei"]  = Ei |> collect
fid["Ej"]  = Ej |> collect
close(fid)

# Parameters to run the second part of the simulation: integrating 
# the optical matrix
kbT   = 0.001
beta  = 1.0/kbT
hw    = 2.0
fermi = 7.0

NE = 2203
NEp = 2401

Gammaf = 1
write  = true # write = true: program writes debug information
distname = "dist.dat"

kappa = 0.0
gph   = 0.6

# Relaxation rates (converted to energy scale)
function Gamma(E,Ep)
    return 0.1
end



Nrows = 12
columns = Array{Float64}(undef, NE, Nrows)
columns[:,1]  = LinRange(Ej[1], Ej[end], NE)
dist, dbd     = sum_dist(Phi,  Ei, Ej, NE, NEp,   fermi,   hw,   beta,   Gammaf,   Gamma, write)
columns[:,2]  = dist
columns[:,3]  = sum_dist(Phi,  Ei, Ej, NE, NEp*3, fermi,   hw,   beta,   Gammaf,   Gamma, false)[2]
columns[:,4]  = sum_dist(Phi,  Ei, Ej, NE, NEp,   fermi+1, hw,   beta,   Gammaf,   Gamma, false)[2]
columns[:,5]  = sum_dist(Phi,  Ei, Ej, NE, NEp,   fermi,   hw+1, beta,   Gammaf,   Gamma, false)[2]
columns[:,6]  = sum_dist(Phi,  Ei, Ej, NE, NEp,   fermi,   hw,   beta*2, Gammaf,   Gamma, false)[2]
columns[:,7]  = sum_dist(Phi,  Ei, Ej, NE, NEp,   fermi,   hw,   beta,   Gammaf*2, Gamma, false)[2]
columns[:,8]  = sum_dist(Phi2, Ei, Ej, NE, NEp,   fermi,   hw,   beta,   Gammaf,   Gamma, false)[2]
columns[:,9]  = sum_dist(Phi3, Ei, Ej, NE, NEp,   fermi,   hw,   beta,   Gammaf,   Gamma, false)[2]

e, Ne, h, dbh = sum_hcg( Phi,  Ei, Ej, NE, NEp,   fermi,   hw,   beta,   kappa, gph, write)
columns[:,10] = Ne
columns[:,11] = sum_hcg(Phi2,  Ei, Ej, NE, NEp,   fermi,   hw,   beta,   kappa, gph, false)[2]
columns[:,12] = sum_hcg(Phi3,  Ei, Ej, NE, NEp,   fermi,   hw,   beta,   kappa, gph, false)[2]

# Write the debug information to file
write_hcg_to_hdf( dbh, "debug_hcg.h5")
write_dist_to_hdf(dbd, "debug_dist.h5")

writedlm(distname, columns)
