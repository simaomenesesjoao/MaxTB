using DelimitedFiles

masterdir = "../"
include(masterdir * "libs/slater_koster.jl")
include(masterdir * "libs/shape_lib.jl")
include(masterdir * "libs/potential.jl")
include(masterdir * "libs/cheb.jl")
include(masterdir * "libs/sumcheb.jl")
include(masterdir * "libs/dist.jl")
include(masterdir * "libs/hcg.jl")
include(masterdir * "libs/dos.jl")


function lowlevel_dist()
    # Script with full control over all parameters of the simulation


    # Properties of the nanoparticle (see documentation for options)
    shape = "cube"
    mater = "gold"
    rad   = 1.5 # [nm] dimension of nanoparticle

    freq  = 2.4 # [eV] frequency of optical potential
    eps_m = 1.0 # relative dielectric constant of the environment

    tempr = 298 # [K] temperature

    outname = "dist.dat"

    # Relaxation in eV
    function Gamma(E,Ep)
        return 0.1
    end

    N = 200 # number of total Chebyshev polynomials

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
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)

    # Get the dielectric constant from frequency and material
    eps = diel(freq)

    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(rad, shape, a0)

    # Use list of atomic positions to determine the Hamiltonian. H is in KPM units
    H = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B)

    # get the potential. Phi is in units of eV
    Phi = potential_sphere(R, eps, eps_m)
    # filename = "pot.dat"
    # Phi = comsol_read(filename)

    # Build the Chebyshev matrix
    mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, Nk)

    # Transform the Chebyshev matrix into the energy-resolved optical matrix
    # Ei,Ej are in eV, opt is dimensionless
    Ei,Ej,opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc1, perc2)

    # Get the distribution
    dist = sum_dist(opt, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, Gammaf, Gamma, write)


    # Write the data
    columns = Array{Float64}(undef, NE, 2)
    columns[:,1] = dist[1] # energy list
    columns[:,2] = dist[2] # population

    # Write the data
    writedlm(outname, columns)
end


lowlevel_dist()