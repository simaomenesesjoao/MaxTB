using Interpolations
using LinearAlgebra

include("aux.jl")


function resum_mu_dd(mumn_orig, Nx, Ny, A, B, perc1, perc2)
    # mu_mn is a matrix of Chebyshev moments, which needs to be resummed
    # Nx,Ny is the number of energies of the first and second indices
    # A,B are the Hamiltonian rescaling factors
    # B - scale
    # A - shift
    #
    # This function resums the matrix with Delta weights in each index
    # Most costly part of this algorithm is the final matrix product, which
    # scales as NTx*NTy*(Nx + Ny)


    # List of energies in the KPM scale (from -1 to 1)
    xlist = LinRange(-0.99, 0.99, Nx)
    ylist = LinRange(-0.99, 0.99, Ny)

    # Convert from KPM to Ha
    Ei = (xlist.*B).+A 
    Ej = (ylist.*B).+A

    # Convert from Ha to eV
    Ei = Ei.*HaeV
    Ej = Ej.*HaeV

    # truncate the Chebyshev matrix
    NTx, NTy = size(mumn_orig)
    NTx  = Int(trunc(perc1/100*NTx))
    NTy  = Int(trunc(perc2/100*NTy))
    mumn = mumn_orig[1:NTx, 1:NTy]

    # Number of Chebyshev polynomials in each index
    # NTx, NTy = size(mumn)
    NTxlist = 0:1:NTx-1|>collect 
    NTylist = 0:1:NTy-1|>collect 

    # apply Jackson kernel to mu matrix (factor 1/2 in first term is included)
    Jx = jackson_coefs(NTx)
    Jy = jackson_coefs(NTy)
    Jxx, Jyy = meshgrid(Jx, Jy)
    mutmn = mumn.*Jxx.*Jyy 

    # Find the Delta coefficients to resum the series
    NTyy, yy = meshgrid(NTylist, ylist)
    xx, NTxx = meshgrid(xlist, NTxlist)

    Dny = cos.(NTyy.*acos.(yy))./sqrt.(1.0.-yy.^2)
    Dxm = cos.(NTxx.*acos.(xx))./sqrt.(1.0.-xx.^2)

    # Resum the Chebyshev series
    phi = (Dxm*mutmn*Dny) / B^2 / HaeV^2

    return Ei, Ej, phi
end 

function test_conv(phi, Ei, Ej)
    # Test if phi is being well calculated by comparing it to its linear interpolation
    Ni = length(Ei)
    Nj = length(Ej)
    println(Ni,Nj)

    # Interpolate only on half of the elements
    phi_itp = LinearInterpolation((Ei[1:2:end], Ej[1:2:end]), phi[1:2:end,1:2:end])

    # Rebuild optical matrix from interpolating function
    phi1 = [phi_itp(Ei[i],Ej[j]) for i=1:Ni-1, j=1:Nj-1]

    # Estimate difference in the matrices
    dif = norm(phi1 - phi[1:end-1,1:end-1])

    return dif
end




