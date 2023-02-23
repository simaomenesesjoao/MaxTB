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

function lowlevel_bs()
    # Script with full control over all parameters of the band structure calculation
    mater = "gold"
    outname = "bands.dat"

    # Select material
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)

    # Calculate the bandstructure along a predefined path
    bands = get_bands(onsite, first_neighbour, second_neighbour)

    # Save to file
    open(outname, "w") do io
        writedlm(io, bands)
    end

end



lowlevel_bs()
