using DelimitedFiles

masterdir = "../src/"
include(masterdir * "libs/slater_koster.jl")
include(masterdir * "libs/shape_lib.jl")
include(masterdir * "libs/potential.jl")
include(masterdir * "libs/cheb.jl")
include(masterdir * "libs/sumcheb.jl")
include(masterdir * "libs/dist.jl")
include(masterdir * "libs/hcg.jl")



function raw()
    # Script with full control over all parameters of the simulation


    # Properties of the nanoparticle (see documentation for options)
    shape = "cube"
    mater = "3DTB"
    rad   = 1.1 # [nm] dimension of nanoparticle

    freq  = 2.0 # [eV] frequency of optical potential
    eps_m = 1.0 # relative dielectric constant of the environment

    tempr = 298 # [K] temperature

    # hot carrier generation parameters
    kappa = 0.0
    gph   = 0.06 # [eV] broadening
    outname = "hcg.dat"


    N = 400 # number of total Chebyshev polynomials

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




    # Select material
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)

    # Get the dielectric constant from frequency and material
    eps = diel(freq)

    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(rad, shape, a0)

    # Use list of atomic positions to determine the Hamiltonian. H is in KPM units
    H, v = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B)

    # get the potential. Phi is in units of eV
    # Phi = potential_sphere(R, eps, eps_m)

    # Build the Chebyshev matrix
    mumn = compute_mumn!(H, v, NNL, NNR, NTx, NTy, Nk)

    # Transform the Chebyshev matrix into the energy-resolved optical matrix
    # Ei,Ej are in eV, opt is dimensionless
    # Ei,Ej,opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc1, perc2)


    # Get the hot carrier generation rate
    # hcg  = sum_hcg( opt, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, kappa, gph, write)


    # columns = hcg_conv(mumn, A, B, NE, fermi, freq, beta, kappa, gph, write)
    print_conv = true
    hcg = hcg_conv(mumn, A, B, NE, fermi_Ha, freq, beta, kappa, gph, print_conv) 



    # Write the data
    # columns = Array{Float64}(undef, NE, 4)
    # columns[:,1] = hcg[1]  # energy list
    # columns[:,2] = hcg[2]  # hot electron
    # columns[:,3] = hcg[3]  # hot hole

    # Write the data
    writedlm(outname, hcg)
end




raw()

