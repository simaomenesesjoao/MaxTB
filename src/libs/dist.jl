using Interpolations
using HDF5
using LinearAlgebra

include("aux.jl")


function build_mask(E, Ep, fermi_eV, beta, hw, Gammaf, Gamma)
    # Prefactor in the distribution calculation which does not depend on the Hamiltonian
    # Gammaf is a factor which multiplies the Gamma function to test for convergence

    # Precalculate the relaxation rates required
    g1 = Gamma(E,Ep)*Gammaf
    g  = Gamma(E,E)*Gammaf

    # Lorentzians
    den_p = (E-Ep+hw)^2 + g1^2
    den_m = (E-Ep-hw)^2 + g1^2
    lorentz = g1/den_p + g1/den_m

    # Difference of Fermi functions
    df = fermi_dirac(E,beta,fermi_eV) - fermi_dirac(Ep,beta,fermi_eV)

    return -4*df/g*lorentz
end


function sum_dist(phi, Ei, Ej, Nx, Ny, fermi_Ha, hw, beta, Gammaf, Gamma, write)
    # Use optical matrix in energy space to build the hot carrier distribution
    # Everything is in units of eV except the Fermi energy [Ha]
    fermi_eV = fermi_Ha*HaeV

    if Nx%2 != 1
        Nx += 1 # must be odd to use in Simpson's integral
    end 

    phi_itp = LinearInterpolation((Ei, Ej), phi)

    # Get limits
    xi = Ei[1]; xf = Ei[end]; yi = Ej[1]; yf = Ej[end];

    xs = LinRange(xi, xf, Nx)
    ys = LinRange(yi, yf, Ny)
    # XX,YY = meshgrid(xs,ys)
    phi1 = phi_itp(xs, ys)
    mask = [build_mask(x,y, fermi_eV, beta, hw, Gammaf, Gamma) for x in xs, y in ys]

    # Complete function to be integrated over
    K = phi1.*mask

              
    # println("size:", size(K))
   

    # Simpson coefficient array
    dy = ys[2] - ys[1]
    simpson = Array{Float64}(undef, Ny)
    simpson[2:2:end] .= 4
    simpson[3:2:end] .= 2
    simpson[1] = 1
    simpson[end] = 1
    simpson .*= dy

    # Integrate first direction using Simpson: 1 4 2 4 2 4 1
    first_integral = Array{Float64}(undef, Nx)
    for j=1:Nx
        first_integral[j] = dot(simpson, K[j,:])
    end

    # Testing 
    db_pack = []
    if write
        db_pack = [Ei,Ej,xs,ys, phi, phi1, mask, K, first_integral]
    end

    return xs, first_integral, db_pack

end

function write_dist_to_hdf(db_pack, name)
    # Write the debug information to an hdf file
    if length(db_pack) == 0
        return 0
    end

    Ei,Ej,xs,ys, phi, phi1, mask, K, first_integral = db_pack

    fid = h5open(name, "w")
    fid["Ei"]  = Ei |> collect
    fid["Ej"]  = Ej |> collect
    fid["xs"]  = xs |> collect
    fid["ys"]  = ys |> collect
    fid["Phi"] = phi1
    fid["mask"] = mask
    fid["K"] = K
    close(fid)
end



function dist_conv(mumn, A, B, NE, fermi_Ha, freq, beta, Gammaf, Gamma, print_conv)
    # Automatically performs a convergence analysis on the dist function
    
    # Percentage of polynomials to keep
    perc100 = 100
    perc70  = 70

    write   = false

    N1, N2 = size(mumn) # number of chebyshev polynomials

    # Number of points to evaluate the optical matrix in
    Nx = N1*2 + 1
    Ny = N2*2 + 1

    # Calculate the optical matrix in energy space
    Ei, Ej, opt1 = resum_mu_dd(mumn, Nx, Ny, A, B, perc100, perc100)

    # Check convergence with the number of energy points
    Ei1  = Ei[1:2:end]
    opt2 = opt1[1:2:end,:]

    NEp  = Nx
    NEp1 = Nx*2

    dist1 = sum_dist(opt1, Ei,  Ej, NE, NEp,  fermi_Ha, freq, beta, Gammaf, Gamma, write)
    dist2 = sum_dist(opt2, Ei1, Ej, NE, NEp,  fermi_Ha, freq, beta, Gammaf, Gamma, write)
    dist3 = sum_dist(opt1, Ei,  Ej, NE, NEp1, fermi_Ha, freq, beta, Gammaf, Gamma, write)

    distance2 = norm(dist1[2] - dist2[2])
    distance3 = norm(dist1[2] - dist3[2])



    # Check convergence with the number of Chebyshev polynomials
    Ei, Ej, opt2 = resum_mu_dd(mumn, Nx, Ny, A, B, perc70,  perc100)
    Ei, Ej, opt3 = resum_mu_dd(mumn, Nx, Ny, A, B, perc100, perc70)
    Ei, Ej, opt4 = resum_mu_dd(mumn, Nx, Ny, A, B, perc70,  perc70)

    dist4 = sum_dist(opt2, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, Gammaf, Gamma, write)
    dist5 = sum_dist(opt3, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, Gammaf, Gamma, write)
    dist6 = sum_dist(opt4, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, Gammaf, Gamma, write)

    distance4 = norm(dist1[2] - dist4[2])
    distance5 = norm(dist1[2] - dist5[2])
    distance6 = norm(dist1[2] - dist6[2])

    # Check convergence with gph and temperature
    dist7 = sum_dist(opt1, Ei,  Ej, NE, NEp,  fermi_Ha, freq, beta, Gammaf*2, Gamma, write)
    dist8 = sum_dist(opt1, Ei,  Ej, NE, NEp,  fermi_Ha, freq, beta/2, Gammaf, Gamma, write)

    distance7 = norm(dist1[2] - dist7[2])
    distance8 = norm(dist1[2] - dist8[2])

    # Check what happens when Fermi and frequency are changed
    dist9  = sum_dist(opt1, Ei,  Ej, NE, NEp,  fermi_Ha+1/HaeV, freq, beta, Gammaf, Gamma, write)
    dist10 = sum_dist(opt1, Ei,  Ej, NE, NEp,  fermi_Ha, freq+1, beta, Gammaf, Gamma, write)

    if print_conv
        println("Double number of energy points in optical matrix ", distance2)
        println("Double number of energy points in energy integration ", distance3)
        println("70% polynomials in first index ", distance4)
        println("70% polynomials in second index ", distance5)
        println("70% polynomials in both indices ", distance6)
        println("Double Gammaf ", distance7)
        println("Double temperature ", distance8)
    end

    # Save all datasets to one array
    columns = Array{Float64}(undef, NE, 11)
    columns[:,1]  = dist1[1]  # Set of output energies
    columns[:,2]  = dist1[2]  # distribution

    columns[:,3]  = dist2[2]  # double energy points in optical matrix
    columns[:,4]  = dist3[2]  # double energy points in integration
    columns[:,5]  = dist4[2]  # 70% polys in first index
    columns[:,6]  = dist5[2]  # 70% polys in second index
    columns[:,7]  = dist6[2]  # 70% polys in both indices
    columns[:,8]  = dist7[2]  # double Gammaf
    columns[:,9]  = dist8[2]  # double temperature
    columns[:,10] = dist9[2]  # add 1 to Fermi
    columns[:,11] = dist10[2] # add 1 to frequency
    return columns

end
