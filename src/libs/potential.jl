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

function linear_pot(R, Norbitals = 9)
    # Builds the diagonal potential from the set of positions, assuming it is
    # proportional to z (exact for dielectric sphere)

    Nat = size(R,1)     # number of atoms
    N   = Nat*Norbitals # size of Hilbert space

    # Put potential in the form of a sparse matrix
    iidx  = zeros(Int64, N)
    jidx  = zeros(Int64, N)
    value = zeros(ComplexF64, N)

    # Iterate over the atoms and orbitals
    for i=1:Nat
        for j=1:Norbitals
            n = (i-1)*Norbitals + j
            iidx[n]  = n
            jidx[n]  = n
            value[n] = R[i,3] # proportional to z
        end 
    end 

    # build potential as sparse matrix
    phi = sparse(iidx, jidx, value)
    return phi
end



function comsol_read(filename;Norbitals=9,Nomean=true,factor=1.0)
    # Read the potential from a COMSOL file
    f=open(filename,"r")

    # Ignore the first 9 lines (header)
    for i=1:9
        😅=readline(f)
    end

    flines=readlines(f)
    φ = Vector{ComplexF64}()
    mean_aux = 0+0im
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
    # ΦT = sparse(iidx, jidx, conj(value))
    return Φ
end



function write_pot(outname, R, phi)
    # Write the potential and positions to a file with COMSOL format
    
    fw = open(outname, "w")

    # Write header
    write(fw, "% Model:              octa2.mph")
    write(fw, "% Version:            COMSOL 5.5.0.359")
    write(fw, "% Date:               Sep 29 2022, 16:37")
    write(fw, "% Dimension:          3")
    write(fw, "% Nodes:              3134")
    write(fw, "% Expressions:        1")
    write(fw, "% Description:        Electric potential")
    write(fw, "% Length unit:        nm")
    write(fw, "x               y               z               V (V)")

    # Iterate over the positions and write to file
    Nat = length(R)   # find number of atoms
    N   = length(phi) # find size of Hilbert space
    if N % Nat != 0
        println("Number of atoms does not divide size of Hilbert space. Exiting")
        exit()
    end
    Norbitals = N ÷ Nat # integer division
    for i in 1:Nat

        x = R[i,1]
        y = R[i,2]
        z = R[i,3]
        j = i*Norbitals

        V = phi[j,j]
        Vr = real(V)
        Vi = imag(V)

        # -9.45953772347708E-10+3.100674740383024E-11i
        write(fw,@sprintf("%f\t\t%f\t\t%f\t\t%+f%+f",x,y,z,Vr,Vi)*"i\n")
    end
    close(fw)

end

function potential_sphere(R, eps, eps_m; Norbitals=9)
    # Get the electric potential inside the sphere

    induced = -3*eps_m/(eps + 2*eps_m)
    phi = linear_pot(R, Norbitals)*induced
    return phi
end


# function write_linear_tofile(infile, outfile, eps_m, eps)
    # Get the electric potential inside the sphere
    # Read positions from file and write potential to another file

    # Fetch the dielectric constant for gold (ε) at this frequency
    # and the dielectric constant of the surrounding environment (ε_m)

    # factor = -3*eps_m/(eps + 2*eps_m)
    # linear_pot(infile, outfile, factor)

# end
