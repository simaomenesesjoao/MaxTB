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
