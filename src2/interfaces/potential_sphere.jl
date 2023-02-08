# Inserts a potential proportional to z at each site position

include("potential.jl")

Nargs = size(ARGS)[1]
if Nargs != 3
    println("Number of arguments (", Nargs, ") is wrong. Usage:\n")
    println("julia potential_sphere.jl [positions] [output filename] [frequency]\n")
    println("Example: julia potential_sphere.jl positions.dat potential.dat 1.8\n")
    println("will calculate the electric potential at positions 'positions.dat' into the file 'potential.dat' at a frequency of 1.8 eV.")
    exit()
end

infile_name  = ARGS[1]
outfile_name = ARGS[2]
freq = parse(Float64,ARGS[3])


compute_potential_sphere_tofile(infile_name, outfile_name, freq)
