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


function lowlevel_hcg()
    # Script with full control over all parameters of the simulation


    # Properties of the nanoparticle (see documentation for options)
    shape = "cube"
    mater = "gold"
    outname = "positions.dat"

    # geometry in nm
    length = 1.0
    height = 1.5
    width  = 2.0

    # Select material
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)

    # Generate list of atomic positions
    Elist, Edict, R = generate_shape_FCC(shape, a0, length, height, width)

    print_positions(R, outname)

end


lowlevel_hcg()
