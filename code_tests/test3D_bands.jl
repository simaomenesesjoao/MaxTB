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

function test(mater)

    # Select material from list
    onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)
    bands = get_bands(onsite, first_neighbour, second_neighbour)

    return bands

end

println(length(ARGS))
if length(ARGS) > 0 && ARGS[1] == "run"

    material = "3DTB"
    bands = test(material)

    open("bands.txt", "w") do io
        writedlm(io, bands)
    end

end
