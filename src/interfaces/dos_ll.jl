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

function lowlevel_dos()
    # Script with full control over all parameters of the simulation

    N  = 300 # number of total Chebyshev polynomials
    Nk = 2   # number of random vectors

    perc = 100   # percentage of Chebyshev moments to keep
    minE_Ha = -1 # smallest energy
    maxE_Ha = +1 # largest energy
    NE = 1000    # number of energies
    flag = 2     # flag=2 means random vectors

    mater = "gold"
    shape = "cube"
    rad = 5.1 # [nm] dimension of nanoparticle
    outname = "dos.dat"

    # Select material
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)


    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(rad, shape, a0)
    println("number of atoms", length(R))

    # Use list of atomic positions to determine the Hamiltonian. H is in KPM units
    H = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B)

    # Build the Chebyshev vector
    mu = doscompute_mu!(H, N, Nk, flag)

    # Resum the Chebyshev vector into the DoS
    dos = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)

    # Save to file
    open(outname, "w") do io
        writedlm(io, dos)
    end



end



lowlevel_dos()
