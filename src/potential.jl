using DataStructures
using Interpolations
using LinearAlgebra
using SparseArrays
using MKLSparse
using Dates
using Random
using MKL
using Printf
using NearestNeighbors

include("dielectric_data.jl")

# KPM shift and scale
# goldA=0.5
# goldB=0.6

# Important constants
HaeV = 27.211386245988       # [eV]        Hartree in eV
boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant


function linear_pot(filename,outname,factor=1)
    # Read the atomic positions and print out a potential which is proportional to z
    # with a proportionality constant of 'factor'. In the case of a dielectric sphere, 
    # this factor can be used as the exact electric potential factor
    
    fr = open(filename,"r")
    fw = open(outname,"w")

    # Copy the first 8 lines (header) to output file
    for i=1:8
        line = readline(fr)
        write(fw,line*"\n")
    end

    # The 9th line gets modified to include "V (V)"
    line = readline(fr)
    write(fw,@sprintf("x\t\ty\t\tz\t\tV (V)\n"))

    # Iterate over the positions and write to file
    flines=readlines(fr)
    for i in flines
        j = split(i)
        x = parse(Float64,j[1])
        y = parse(Float64,j[2])
        z = parse(Float64,j[3])
        V = z*factor
        Vr = real(V)
        Vi = imag(V)
        # -9.45953772347708E-10+3.100674740383024E-11i
        write(fw,@sprintf("%f\t\t%f\t\t%f\t\t%+f%+f",x,y,z,Vr,Vi)*"i\n")
    end
    close(fr)
    close(fw)

end



function compute_potential_sphere_tofile(infile, outfile, freq)
    # Get the electric potential inside the sphere
    # Read positions from file and write potential to another file

    # Fetch the dielectric constant for gold (ε) at this frequency
    # and the dielectric constant of the surrounding environment (ε_m)
    ε = get_dielectric(freq, "gold")
    ε_m = eps_m

    a_0 = 1

    # factor m
    factor = -3*ε_m/(ε + 2*ε_m)*a_0/2
    linear_pot(infile, outfile, factor)

end


function compute_potential_sphere(R,Edict)
    # Get the electric potential inside the sphere
    #
    a,b=size(R)
    value_aux=zeros(ComplexF64,a)
    # Potential is proportional to z
    for i=1:a
        RR=R[i,:]
        value_aux[i]=-3*eps_m/(eps_pd+2*eps_m)*RR[3]*a_0/2
    end 

    iidx=1:length(Edict)*9|>collect
    jidx=1:length(Edict)*9|>collect

    value=zeros(ComplexF64,length(Edict)*9)
    for i=1:length(Edict)
        for j=1:9
            value[(i-1)*9+j]=value_aux[i]
        end 
    end 
    Phi=sparse(iidx,jidx,value)
    return Phi
end


function comsol_read(filename;Norbitals=9,Nomean=true,factor=1.0)
    # Read the potential from a COMSOL file
    f=open(filename,"r")

    # Ignore the first 9 lines (header)
    for i=1:9
        😅=readline(f)
    end

    flines=readlines(f)
    φ=Vector{ComplexF64}()
    mean_aux=0+0im
    for i in flines
        if length(i)>3
        j=split(i)
        φi=parse(ComplexF64,j[4])*factor
        push!(φ,φi)
        mean_aux=mean_aux+φi
        end
    end

    # Remove mean, reduces fluctuations
    if Nomean
        mean=mean_aux/length(φ)
        φ=φ.-mean
    end 

    # Put potential in the form of a sparse matrix
    iidx=zeros(Int64,length(φ)*Norbitals)
    jidx=zeros(Int64,length(φ)*Norbitals)
    value=zeros(ComplexF64,length(φ)*Norbitals)
    for i=1:lastindex(φ)
        for j=1:Norbitals
            iidx[(i-1)*Norbitals+j]=(i-1)*Norbitals+j
            jidx[(i-1)*Norbitals+j]=(i-1)*Norbitals+j
            value[(i-1)*Norbitals+j]=φ[i]
        end 
    end 

    Φ  = sparse(iidx, jidx, value)
    ΦT = sparse(iidx, jidx, conj(value))
    return Φ
end

