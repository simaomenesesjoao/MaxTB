# Inserts a potential proportional to z at each site position

include("potential.jl")

Nargs = size(ARGS)[1]
if Nargs != 3
    println("Number of arguments (", Nargs, ") is wrong. Usage:\n")
    println("julia potential_sphere.jl [positions file] [potential filename] [frequency]\n")
    println("Example: julia potential_sphere.jl positions.dat potential.dat 2.0\n")
    println("will calculate the electrostatic potential for every point inside the file positions.dat, and will output a very similar file but with a fourth column containing the complex potential")
    exit()
end

infile_name  = ARGS[1]
outfile_name = ARGS[2]
freq = parse(Float64,ARGS[3])

compute_potential_sphere_tofile(infile_name, outfile_name, freq)
println("Finished computing potential.")
