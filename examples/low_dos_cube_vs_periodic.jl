using DelimitedFiles

masterdir = "../src/"
include(masterdir * "libs/slater_koster.jl")
include(masterdir * "libs/shape_lib.jl")
include(masterdir * "libs/potential.jl")
include(masterdir * "libs/cheb.jl")
include(masterdir * "libs/sumcheb.jl")
include(masterdir * "libs/dist.jl")
include(masterdir * "libs/hcg.jl")
include(masterdir * "libs/dos.jl")

function dos()
    N  = 300 # number of total Chebyshev polynomials
    Nk = 2   # number of random vectors

    perc = 100   # percentage of Chebyshev moments to keep
    minE_Ha = -1 # smallest energy
    maxE_Ha = +1 # largest energy
    NE = 1000    # number of energies
    flag = 2     # flag=2 means random vectors

    # Select material
    mater = "gold"
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)


    shape = "periodic"
    L = 8 # number of crystal planes for the PBC
    rad = (L-1+0.01)*a0/2 # [nm] dimension of nanoparticle

    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(rad, shape, a0)
    println("number of atoms", length(R))

    # Use list of atomic positions to determine the Hamiltonian. H is in KPM units
    H = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B, L1=L, L2=L, L3=L, periodic=true)

    # Build the Chebyshev vector
    mu = doscompute_mu!(H, N, Nk, flag)

    # Resum the Chebyshev vector into the DoS
    dos_periodic = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)





    shape = "cube"
    rad = 4.1 # size of the cube

    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(rad, shape, a0)
    println("number of atoms", length(R))

    # Use list of atomic positions to determine the Hamiltonian. H is in KPM units
    H = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B)

    # Build the Chebyshev vector
    mu = doscompute_mu!(H, N, Nk, flag)

    # Resum the Chebyshev vector into the DoS
    dos_cube = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)

    open("dos_periodic.dat", "w") do io
        writedlm(io, dos_periodic)
    end
    open("dos_cube.dat", "w") do io
        writedlm(io, dos_cube)
    end



end




dos()

