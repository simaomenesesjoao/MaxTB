# Include the materials
include("../materials/gold.jl")
include("../materials/silver.jl")
include("../materials/al.jl")
include("../materials/palladium.jl")
include("../materials/cube3D.jl")

function tightbinding(material)
    if material == "gold"
        return gold()
    elseif material == "silver"
        return silver()
    elseif material == "palladium"
        return palladium()
    elseif material == "aluminium"
        return aluminium()
    elseif material == "3DTB"
        return cube()
    else
        println("Material not supported")
        return -1
    end
end
