using Interpolations
using HDF5
using LinearAlgebra



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
