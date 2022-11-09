shape_name = ARGS[1]
Rmax       = parse(Float64, ARGS[2])
rescale    = parse(Bool,    ARGS[3])

shapes = ["octahedron", "cube", "rhombic", "sphere"]

if !(shape_name in shapes)
    println("Shape "*shape_name*" not supported.")
    exit()
end

include("shape_lib.jl")

# Generate the positions of the atoms inside the nanoparticle in nanometers
# Rmax is the length of the solid's edge in nanometers
l = Rmax/(a_0/2) # number of atomic (100) planes
Elist, Edict, R = generate_shape(l, shape_name)
println("Finished creating the atomic positions. Number of atoms: ", length(Elist))

Rstr = string(Rmax)
a = split(Rstr, ".")
Rstr = a[1]*"p"*a[2]
# filename = "positions_"*shape_name*"_rad"*Rstr*".dat"
# filename = shape_name*"/positions_rad"*Rstr*".dat"
filename = "positions.dat"
println("filename: ", filename)
# print_positions(R, filename)
print_positions_comsol(R, filename, rescale=rescale)
