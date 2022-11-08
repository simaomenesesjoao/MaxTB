using HDF5
using DelimitedFiles
using Printf

goldA=0.5
goldB=0.6

include("aux.jl")

# Important constants
HaeV = 27.211386245988       # [eV]        Hartree in eV
boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant
#a_0 = 7.291                  # Lattice constant

# Important parameters
#Fermi = 0.2595               # [Ha]     Fermi energy
T     = 298                  # [Kelvin] room temperature 
beta  = 1/(boltzmann*T/HaeV) # [Ha⁻¹]   inverse temperature

# Modification to the spectral lines
kappa = 0.0
gph   = 0.06/HaeV            # [Ha] broadening


Nargs   = size(ARGS)[1]
if Nargs != 6
    println("Number of arguments is wrong. h5name, outname, freq, Fermi, NE, %")
end

h5name  = ARGS[1] # name of the file with the Chebyshev moments
outname = ARGS[2] # name of the output file with the hot carrier generation rate
hw      = parse(Float64, ARGS[3]) # frequency (eV)
Fermi   = parse(Float64, ARGS[4]) # Fermi energy (Hartree)
N1      = parse(Int,     ARGS[5]) # number of energies for the energy integration
perc    = parse(Float64, ARGS[6]) # percentage of the Chebyshev polynomials to use

println("Arguments read from command line: ", h5name, outname, hw, Fermi, N1, perc)

# Process the input
N2 = N1
hw = hw/HaeV # convert eV to Hartree

# read from the hdf5 file
fid = h5open(h5name, "r")
mumn = read(fid["moments"])
close(fid)

# truncate the Chebyshev matrix
n,m = size(mumn)
n = Int(trunc(perc/100*n))
m = Int(trunc(perc/100*m))
mumn = mumn[1:n, 1:m]

println("Resumming Chebyshev moments")
Ej,Nelist,Nhlist=compute_eh_list(mumn,N1,N2,hw,goldA,goldB,Fermi,beta,kappa,gph)

f3=open(outname,"w")
for i=1:lastindex(Ej)
    write(f3,@sprintf(" %+1.15E",Ej[i]))
    write(f3,@sprintf(" %+1.15E",Nelist[i]))
    write(f3,@sprintf(" %+1.15E",Nhlist[i]))
    write(f3,"\n")
end


