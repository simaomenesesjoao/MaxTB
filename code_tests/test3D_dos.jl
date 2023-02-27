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

function test(L,N, mater)
    # L: number of crystal planes for the PBC
    # N: number of total Chebyshev polynomials
    # mater: material

    # Properties of the nanoparticle (see documentation for options)
    shape = "periodic"


    # Select material from list
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)

    rad = (L-1+0.01)*a0/2 # [nm] dimension of nanoparticle

    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(rad, shape, a0)
    println("number of atoms: ", length(R))

    # Use list of atomic positions to determine the Hamiltonian. H is in KPM units
    H,v = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B, L1=L, L2=L, L3=L, periodic=true)

    # Build the Chebyshev matrix
    flag = 2    # use random vectors for the calculation
    Nk   = 10   # number of random vectors
    mu   = doscompute_mu!(H, N, Nk, flag)

    # Get the DoS from the Chebyshev matrix
    perc = 100 # percentage of polynomials to keep
    minE_Ha = -1 # smallest energy
    maxE_Ha = +1 # largest energy 
    NE = 1000 # number of energy points
    dos = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)


    return dos




end

println(length(ARGS))
if length(ARGS) > 0 && ARGS[1] == "run"
    L = 8
    N = 200
    mater = "3DTB"
    dos = test(L, N, mater)

    open("dos.txt", "w") do io
        writedlm(io, transpose(dos))
    end

end
