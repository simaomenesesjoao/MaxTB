shape_name = ARGS[1]
shapes = ["octahedron", "cube", "rhombic"]

if !(shape_name in shapes)
    println("Shape "*shape_name*" not supported.")
    exit()
end

Rmax = parse(Float64, ARGS[2])

include("shape_lib.jl")

# Generate the positions of the atoms inside the nanoparticle
Elist, Edict, R = generate_shape(Rmax, shape_name)
println("Finished creating the atomic positions. Number of atoms: ", length(Elist))

filename = shape_name*"/positions_rad"*string(Rmax)*".dat"
# print_positions(R, filename)
print_positions_comsol(R, filename)
