# This is an interface which builds the optical matrix in energy 
# space and integrates it to produce the hcg and the distribution, 
# but also produces some debugging information
# The algorithms scale as M²(N1 + N2) + NEp*NE
using HDF5
using DelimitedFiles
using Printf

include("sumcheb.jl")
include("dist.jl")
include("hcg.jl")

HaeV = 27.211386245988       # [eV]        Hartree in eV
boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant

Nargs = size(ARGS)[1]
if Nargs != 8
    println("Number of arguments is wrong. Should be 8")
end

h5name  = ARGS[1]
outname = ARGS[2] # name of the output file with the electron distribution
hw      = parse(Float64, ARGS[3])
fermi   = parse(Float64, ARGS[4])*HaeV
N1      = parse(Int,     ARGS[5])
N2      = parse(Int,     ARGS[6])
NE      = parse(Int,     ARGS[7])
NEp     = parse(Int,     ARGS[8])

tempr   = 298
kbT     = tempr*boltzmann
beta    = 1.0/kbT

fid = h5open(h5name, "r")
mumn = read(fid["moments"])
close(fid)

A = 0.5
B = 1.0
Gammaf = 1.0
kappa = 0.0
gph   = 0.06/HaeV            # [Ha] broadening

perc1 = 100
perc2 = 70
perc3 = 40
Eia, Eja, Phi1a = resum_mu_dd(mumn, N1*2, N2*2, A, B, perc1, perc1)
Ei,  Ej,  Phi1  = resum_mu_dd(mumn, N1,   N2,   A, B, perc1, perc1)
Ei,  Ej,  Phi2  = resum_mu_dd(mumn, N1,   N2,   A, B, perc2, perc2)
Ei,  Ej,  Phi3  = resum_mu_dd(mumn, N1,   N2,   A, B, perc3, perc3)

Ei  *= HaeV
Ej  *= HaeV
Eia *= HaeV
Eja *= HaeV

write = false
Nrows = 16
columns = Array{Float64}(undef, NE, Nrows)
columns[:,1]  = LinRange(Ej[1], Ej[end], NE)

function Gamma(E,Ep)
    return 0.1
end


# Distribution (change the optical matrix)
columns[:,2]  = sum_dist(Phi1a, Eia, Eja, NE, NEp,   fermi,   hw,   beta,   Gammaf,   Gamma, write)[2]
columns[:,3]  = sum_dist(Phi1,  Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   Gammaf,   Gamma, write)[2]
columns[:,4]  = sum_dist(Phi2,  Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   Gammaf,   Gamma, write)[2]
columns[:,5]  = sum_dist(Phi3,  Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   Gammaf,   Gamma, write)[2]

# Distribution (fix optical matrix, change all other parameters)
columns[:,6]  = sum_dist(Phi1,  Ei,  Ej,  NE, NEp*3, fermi,   hw,   beta,   Gammaf,   Gamma, write)[2]
columns[:,7]  = sum_dist(Phi1,  Ei,  Ej,  NE, NEp,   fermi+1, hw,   beta,   Gammaf,   Gamma, write)[2]
columns[:,8]  = sum_dist(Phi1,  Ei,  Ej,  NE, NEp,   fermi,   hw+1, beta,   Gammaf,   Gamma, write)[2]
columns[:,9]  = sum_dist(Phi1,  Ei,  Ej,  NE, NEp,   fermi,   hw,   beta*2, Gammaf,   Gamma, write)[2]
columns[:,10] = sum_dist(Phi1,  Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   Gammaf*2, Gamma, write)[2]

# Hot carrier generation (change the optical matrix)
columns[:,11] = sum_hcg(Phi1,   Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   kappa,    gph,   write)[2]
columns[:,12] = sum_hcg(Phi2,   Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   kappa,    gph,   write)[2]
columns[:,13] = sum_hcg(Phi3,   Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   kappa,    gph,   write)[2]

# Hot carrier generation (fix optical matrix, change all other parameters)
columns[:,14] = sum_hcg(Phi1,   Ei,  Ej,  NE, NEp,   fermi,   hw,   beta,   kappa*2,  gph*2, write)[2]
columns[:,15] = sum_hcg(Phi1,   Ei,  Ej,  NE, NEp,   fermi,   hw,   beta*2, kappa,    gph,   write)[2]
columns[:,16] = sum_hcg(Phi1,   Ei,  Ej,  NE, NEp*2, fermi,   hw,   beta,   kappa,    gph,   write)[2]

# Write the data
writedlm(outname, columns)
