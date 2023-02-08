filename   = ARGS[1]
shape_name = ARGS[2]
Rmax       = parse(Float64, ARGS[3])
rescale    = parse(Bool,    ARGS[4])


Nargs = size(ARGS)[1]
if Nargs != 4
    println("Number of arguments (", Nargs, ") is wrong. Usage:\n")
    println("julia generate_positions.jl [filename] [shape name] [length] [rescale].\n")
    println("Example: julia generate_positions.jl positions.dat sphere 2.0 false\n")
    println("will output the atomic positions inside the nanoparticle specified by the shape (sphere) and length (2.0). In the case of a sphere, the length represents the radius. Other shapes: octahedron, cube and rhombic (rhombic dodecahedron). In this case, the length represents the edge length of these solids.")
    exit()
end

shapes = ["octahedron", "cube", "rhombic", "sphere"]
if !(shape_name in shapes)
    println("Shape "*shape_name*" not supported. Use sphere, octahedron, cube or rhombic")
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
# filename = "positions.dat"
println("filename: ", filename)
# print_positions(R, filename)
print_positions_comsol(R, filename, rescale=rescale)
