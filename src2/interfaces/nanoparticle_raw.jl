using DelimitedFiles

include("slater_koster.jl")
include("shape_lib.jl")
include("potential.jl")
include("cheb.jl")
include("sumcheb.jl")
include("dist.jl")
include("hcg.jl")

boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant
HaeV = 27.211386245988       # [eV]        Hartree in eV


function raw()
    # Script with full control over all parameters of the simulation


    # Properties of the nanoparticle
    shape = "cube"
    mater = "gold"
    rad   = 1.5

    # Dielectric properties
    freq  = 2.4 # frequency of optical potential
    eps_m = 1.0

    N = 200 # number of total Chebyshev polynomials

    # Temperature
    tempr   = 298

    # hot carrier generation parameters
    kappa = 0.0
    gph   = 0.06/HaeV            # [Ha] broadening
    outname = "hcg.dat"

    # Relaxation
    function Gamma(E,Ep)
        return 0.1
    end




    # number of Chebyshev polynomials per block
    NTx = 10
    NTy = 10

    # number of blocks
    NNL = N÷NTx + 1
    NNR = N÷NTy + 1
    Nk  = 5 # number of random vectors


    # Number of integration points in the energy discretization
    Nx = N*2
    Ny = Nx  

    # Percentage of Chebyshev polynomials to keep
    perc1 = 100
    perc2 = perc1

    NE  = N*2+1 # Number of energy points in the integration
    NEp = N*2+1 # Number of output energy points

    # Temperature
    kbT  = tempr*boltzmann # eV
    beta = 1.0/kbT

    write  = false # write debug information
    Gammaf = 1 # prefactor to the relaxation function 




    # Select material
    onsite, first_neighbour, second_neighbour, A, B, fermi, a0, diel = tightbinding(mater)

    # Get the dielectric constant from frequency and material
    eps = diel(freq)

    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(rad, shape, a0)

    # Use list of atomic positions to determine the Hamiltonian
    H = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B)

    # get the potential
    Phi = potential_sphere(R, eps, eps_m)

    # Build the Chebyshev matrix
    mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, Nk)

    # Transform the Chebyshev matrix into the energy-resolved optical matrix
    Ei,Ej,opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc1, perc2)
    Ei *= HaeV # convert to eV
    Ej *= HaeV # convert to eV

    # Get the distribution
    dist = sum_dist(opt, Ei, Ej, NE, NEp, fermi, freq, beta, Gammaf, Gamma, write)

    # Get the hot carrier generation rate
    hcg  = sum_hcg( opt, Ei, Ej, NE, NEp, fermi, freq, beta, kappa, gph, write)


    # columns = hcg_conv(mumn, A, B, NE, fermi, freq, beta, kappa, gph, write)
    print_conv = true
    hcg_conv(mumn, A, B, NE, fermi, freq, beta, kappa, gph, print_conv) 



    # Write the data

    columns = Array{Float64}(undef, NE, 4)
    columns[:,1] = hcg[1]  # energy list
    columns[:,2] = hcg[2]  # hot electron
    columns[:,3] = hcg[3]  # hot hole
    columns[:,4] = dist[2] # population

    # Write the data
    writedlm(outname, columns)
end




raw()

