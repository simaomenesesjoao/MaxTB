include("aux.jl")

function build_mask_hcg(E, Ep, fermi, beta, hw, kappa, gph)
    # Prefactors in the hot carrier calculation which does not depend on the Hamiltonian
    Nx = length(E)
    Ny = length(Ep)

    function getgi(Ei,gph,kappa,fermi)
        ans=(kappa.*(Ei.-fermi).^2).+gph
        return ans 
    end 

    gi = getgi(E,gph,kappa,fermi)
    gj = getgi(Ep,gph,kappa,fermi)

    # Integrating step: kernel
    factor = zeros(Float64,Nx,Ny)
    hole_factor = zeros(Float64,Nx,Ny)
    for i=1:Nx 
        Ei = E[i]
        fi = 1.0/(1.0 + exp(beta*(Ei-fermi)))

        for j=1:Ny 
            Ej = Ep[j]
            fj = 1.0/(1.0 + exp(beta*(Ej-fermi)))

            gij  = gi[i] + gj[j] 
            amp  = 1.0/sqrt(2.0*pi*gij^2) # amplitude of the gaussian
            d_ij = amp*exp(-(hw+(Ej-Ei ) )^2/2.0/gij^2)
            d_ji = amp*exp(-(hw-(Ej-Ei ) )^2/2.0/gij^2)
            
            factor[i,j]      = d_ji*fi*(1.0-fj)
            hole_factor[i,j] = d_ij*fj*(1.0-fi)
        end 
    end 

    return factor, hole_factor

end


function sum_hcg(phi,Ei,Ej,Nx,Ny,fermi, hw, beta, kappa, gph, debug)
    # Use optical matrix in energy space to build the hot carrier 
    # generation rate. Everything is in units of eV

    Nx,Ny = Ny,Nx
    # Build interpolating function
    phi_itp = LinearInterpolation((Ei, Ej), phi)

    # Get new arrays using interpolation
    xi = Ei[1]; xf = Ei[end]; yi = Ej[1]; yf = Ej[end]
    xs = LinRange(xi, xf, Nx)
    ys = LinRange(yi, yf, Ny)

    phi1 = phi_itp(xs, ys)
    mask_el, mask_ho = build_mask_hcg(xs,ys, fermi, beta, hw, kappa, gph)

    # Apply the masks to the optical matrices in energy space
    Ke = phi1.*mask_el
    Kh = phi1.*mask_ho
    
    Nelist = zeros(Ny)
    Nhlist = zeros(Ny)
    
    # Integration step
    dx = xs[2] - xs[1]
    for i=1:Ny 
        for j=1:Nx
            Nelist[i] = Nelist[i] + Ke[j,i]*dx
            Nhlist[i] = Nhlist[i] + Kh[j,i]*dx
        end 
    end 

    factor = 2*pi*4/pi^2
    Nelist = Nelist.*factor
    Nhlist = Nhlist.*factor

    # debugging option. Save all the energy matrices to an hdf file
    db_pack = []
    if debug == 2
        db_pack = [Ei,Ej,xs,ys, phi, phi1, mask_el, mask_ho, Ke, Kh, Nelist, Nhlist]
    end
    
    return ys, Nelist, Nhlist, db_pack
end 



function write_hcg_to_hdf(db_pack, name)
    # Write the debug information to an hdf file
    if length(db_pack) == 0
        return 0
    end
    Ei,Ej,xs,ys, phi, phi1, mask_el, mask_ho, Ke, Kh, Nelist, Nhlist = db_pack

    fid = h5open(name, "w")

    fid["Phi"] = phi
    fid["Ei"]  = Ei |> collect
    fid["Ej"]  = Ej |> collect

    fid["Phi_itp"] = phi1
    fid["mask_el"] = mask_el
    fid["mask_ho"] = mask_ho
    fid["Ke"] = Ke
    fid["Kh"] = Kh
    fid["xs"]  = xs |> collect
    fid["ys"]  = ys |> collect

    fid["Nelist"] = Nelist
    fid["Nhlist"] = Nhlist

    close(fid)
end




function hcg_auto(mumn, A, B, NE, fermi, freq, beta, kappa, gph, write)
    # Calculate the hot carrier generation rate with sensible default parameters

    perc100 = 100
    write   = false

    N1, N2 = size(mumn) # number of chebyshev polynomials

    # Number of points to evaluate the optical matrix in
    Nx  = N1*2 + 1
    Ny  = N2*2 + 1

    NEp = Nx # Number of integration energies

    # Calculate the optical matrix in energy space
    Ei, Ej, opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc100, perc100)

    # Use the optical matrix to obtain the hot carrier generation rate
    hcg = sum_hcg(opt, Ei,  Ej, NE, NEp,  fermi, freq, beta, kappa, gph, write)

    return hcg
end


function hcg_conv(mumn, A, B, NE, fermi, freq, beta, kappa, gph, print_conv)
    # Automatically performs a convergence analysis on the hcg function
    
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

    hcg1 = sum_hcg(opt1, Ei,  Ej, NE, NEp,  fermi, freq, beta, kappa, gph, write)
    hcg2 = sum_hcg(opt2, Ei1, Ej, NE, NEp,  fermi, freq, beta, kappa, gph, write)
    hcg3 = sum_hcg(opt1, Ei,  Ej, NE, NEp1, fermi, freq, beta, kappa, gph, write)

    dist2 = norm(hcg1[2] - hcg2[2])
    dist3 = norm(hcg1[2] - hcg3[2])



    # Check convergence with the number of Chebyshev polynomials
    Ei, Ej, opt2 = resum_mu_dd(mumn, Nx, Ny, A, B, perc70,  perc100)
    Ei, Ej, opt3 = resum_mu_dd(mumn, Nx, Ny, A, B, perc100, perc70)
    Ei, Ej, opt4 = resum_mu_dd(mumn, Nx, Ny, A, B, perc70,  perc70)

    hcg4 = sum_hcg(opt2, Ei, Ej, NE, NEp, fermi, freq, beta, kappa, gph, write)
    hcg5 = sum_hcg(opt3, Ei, Ej, NE, NEp, fermi, freq, beta, kappa, gph, write)
    hcg6 = sum_hcg(opt4, Ei, Ej, NE, NEp, fermi, freq, beta, kappa, gph, write)

    dist4 = norm(hcg1[2] - hcg4[2])
    dist5 = norm(hcg1[2] - hcg5[2])
    dist6 = norm(hcg1[2] - hcg6[2])

    # Check convergence with gph and temperature
    hcg7 = sum_hcg(opt1, Ei,  Ej, NE, NEp,  fermi, freq, beta, kappa, gph*2, write)
    hcg8 = sum_hcg(opt1, Ei,  Ej, NE, NEp,  fermi, freq, beta/2, kappa, gph, write)

    dist7 = norm(hcg1[2] - hcg7[2])
    dist8 = norm(hcg1[2] - hcg8[2])

    # Check what happens when Fermi and frequency are changed
    hcg9  = sum_hcg(opt1, Ei,  Ej, NE, NEp,  fermi+1, freq, beta, kappa, gph, write)
    hcg10 = sum_hcg(opt1, Ei,  Ej, NE, NEp,  fermi, freq+1, beta, kappa, gph, write)

    if print_conv
        println("Double number of energy points in optical matrix", dist2)
        println("Double number of energy points in energy integration", dist3)
        println("70% polynomials in first index", dist4)
        println("70% polynomials in second index", dist5)
        println("70% polynomials in both indices", dist6)
        println("Double gph", dist7)
        println("Double temperature", dist8)
    end

    # Save all datasets to one array
    columns = Array{Float64}(undef, NE, 12)
    columns[:,1] = hcg1[1] # Set of output energies
    columns[:,2] = hcg1[2] # Hot electrons
    columns[:,3] = hcg1[3] # Hot holes

    columns[:,4]  = hcg2[2]  # double energy points in optical matrix
    columns[:,5]  = hcg3[2]  # double energy points in integration
    columns[:,6]  = hcg4[2]  # 70% polys in first index
    columns[:,7]  = hcg5[2]  # 70% polys in second index
    columns[:,8]  = hcg6[2]  # 70% polys in both indices
    columns[:,9]  = hcg7[2]  # double gph
    columns[:,10] = hcg8[2]  # double temperature
    columns[:,11] = hcg9[2]  # add 1 to Fermi
    columns[:,12] = hcg10[2] # add 1 to frequency
    return columns

end
