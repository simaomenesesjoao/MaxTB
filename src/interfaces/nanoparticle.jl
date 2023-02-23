include("../libs/slater_koster.jl")
include("../libs/shape_lib.jl")
include("../libs/potential.jl")
include("../libs/cheb.jl")
include("../libs/sumcheb.jl")
include("../libs/dist.jl")
include("../libs/hcg.jl")
include("../libs/dos.jl")

# Include the gold and palladium datasets
# include("../materials/gold.jl")
# include("../materials/palladium.jl")
# include("../materials/cube3D.jl")

# function tightbinding(material)
    # if material == "gold"
        # return gold()
    # elseif material == "palladium"
        # return palladium()
    # elseif material == "cube"
        # return cube()
    # else
        # println("Material not supported")
        # return -1
    # end
# end

using DelimitedFiles

function hamiltonian_builder_predefined(shape::AbstractString, rad_nm::Float64, material::AbstractString)
    # Select from a predefined selection of materials and shapes. This function checks the input for errors
    
    # Check if the material selected is implemented
    materials_implemented = ["gold", "palladium"]
    if !(material in materials_implemented)
        println("Material ", material, " not implemented yet. Exiting.\n")
        exit()
    end

    # Check if the shape has been implemented
    shapes_implemented = ["sphere", "octahedron", "cube", "dodecahedron"]
    if !(shape in shapes_implemented)
        println("Shape ", shape, " not implemented yet. Exiting.\n")
        exit()
    end

    if rad_nm <= 0
        println("The radius has to be positive. Exiting.\n")
    end

    description = "predefined"

    return [description, shape, rad_nm, material]
end




function potential_builder_sphere(freq_eV::Float64; eps = "default", eps_m::Float64 = 1.0)

    if freq_eV <= 0.0
        println("The frequency has to be positive. Exiting.\n")
        exit()
    end

    description = "sphere"
    return [description, freq_eV, eps, eps_m]
end

function potential_builder_import(freq_eV::Float64, filename::AbstractString)
    if freq_eV <= 0.0
        println("The frequency has to be positive. Exiting.\n")
        exit()
    end
    description = "import"
    return [description, freq_eV, filename]

end

function calculation_builder_dos(;NE::Int64=1000, minE_Ha::Float64=-2, maxE_Ha::Float64=2, N::Int64=200, NK::Int=5, name::AbstractString="dos.dat")
    description = "dos"
    return [description, NE, minE_Ha, maxE_Ha, N, NK, name]
end

function calculation_builder_bands(;name::AbstractString="bands.dat")
    description = "band structure"
    return [description, name]
end


function calculation_builder_hcg(;temp::Float64=298.0, gph::Float64=0.01, N::Int64=200, NR::Int64=5, NB::Int64=10, NE::Int64=1000, conv::Bool=false, name::AbstractString="hcg.dat")

    # Check if the quantities make sense
    quantities = [temp, gph, N, NR, NB, NE]
    varnames = ["Temperature", "Lifetime", "Number of Chebyshev polynomials", "Number of random vectors", "Number of blocks", "Number of output energies"]

    for i in 1:6
        q = quantities[i]
        v = varnames[i]
        if q<0
            println(v, " cannot be negative. Exiting.\n")
            exit()
        end
    end

    description = "hcg"
    return [description, temp, gph, N, NR, NB, NE, conv, name]
end



function default(E,Ep)
    return 1.0
end

function calculation_builder_dist(;temp::Float64=298, Gamma=default, N::Int64=200, NR::Int64=5, NB::Int64=10, NE::Int64=1000, conv::Bool=false, name::AbstractString="dist.dat")

    # Check if the quantities make sense
    quantities = [temp, N, NR, NB, NE]
    varnames = ["Temperature", "Number of Chebyshev polynomials", "Number of random vectors", "Number of blocks", "Number of output energies"]

    for i in 1:5
        q = quantities[i]
        v = varnames[i]
        if q<0
            println(v, " cannot be negative. Exiting.\n")
            exit()
        end
    end

    description = "dist"
    return [description, temp, Gamma, N, NR, NB, NE, conv, name]
end











function nanoparticle_run(;ham_builder=[], pot_builder=[], hcg_builder=[], dist_builder=[], dos_builder=[], band_builder=[])
    # Perform the complete workflow to obtain the hot carrier generation/distribution 


    # Check which builders have been specified
    ham_b = length(ham_builder ) > 0
    pot_b = length(pot_builder ) > 0
    hcg_b = length(hcg_builder ) > 0
    dis_b = length(dist_builder) > 0
    dos_b = length(dos_builder ) > 0
    ban_b = length(band_builder) > 0


    need_mat = hcg_b || dis_b || dos_b || ban_b
    need_at  = hcg_b || dis_b || dos_b 
    need_H   = hcg_b || dis_b || dos_b 
    need_V   = hcg_b || dis_b

    if need_H && !ham_b
        println("Hamiltonian required for computation but Hamiltonian builder was not specified. Exiting.\n")
        exit()
    end

    if need_V && !pot_b
        println("Potential required for computation but potential builder was not specified. Exiting.\n")
        exit()
    end


    if need_mat
        mater  = ham_builder[4] 

        # Fetch the tight-binding model parameters
        onsite, first_neighbour, second_neighbour, A, B, fermi_Ha, a0, diel = tightbinding(mater)
    end

    if need_at
        shape  = ham_builder[2]
        rad_nm = ham_builder[3]

        # Generate list of atomic positions
        Elist, Edict, R = generate_shape_FCC(rad_nm, shape, a0)
    end

    # Get Hamiltonian from the predefined set
    if need_H
        ham_desc = ham_builder[1]

        if ham_desc == "predefined"

            # Use list of atomic positions to determine the Hamiltonian
            H = slater_koster_FCC(Elist, Edict, onsite, first_neighbour, second_neighbour, A, B)


        # Get Hamiltonian sparse matrix from file
        elseif ham_desc == "custom"
            ham_fname = ham_builder[2]
            # H, A, B, fermi = get_hami(ham_fname)

        else 
            println("Requested Hamiltonian builder is not supported. Exiting.\n")
            exit()
        end

    end



    if need_V
        # Descriptoin and frequency of the optical potential. Always required
        pot_desc = pot_builder[1]
        freq     = pot_builder[2]

        # Dielectric sphere: pot_builder = [pot_desc, freq, eps, eps_m]
        if pot_desc == "sphere"
            eps = pot_builder[3]
            if eps == "default"
                eps = get_dielectric(freq, mater)
            end
            eps_m = pot_builder[4] # dielectric of the environment
            Phi   = potential_sphere(R, eps, eps_m)

        # Reading from file: pot_builder = [pot_desc, freq, filename]
        elseif pot_desc == "import"
            pot_fname = pot_builder[3]
            Phi = comsol_read(pot_fname)

        else
            print("Requested Potential builder is not supported. Exiting.\n")
            exit()
        end
    end





    if hcg_b
        println("Computing the hot carrier generation rate.")
        hcg_desc = hcg_builder[1] 
        # number of total Chebyshev polynomials
        # description, temperature, lifetime, number of chebyshev polynomials,
        # number of random vectors, number of blocks, convergence, name
        desc, tempr, gph, N, NR, NB, NE, conv, name  = hcg_builder

        # number of Chebyshev polynomials per block
        NTx = NB
        NTy = NB

        # number of blocks
        NNL = N÷NTx + 1
        NNR = N÷NTy + 1

        # Build the Chebyshev matrix
        mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, NR)

        # Get inverse temperature
        kbT     = tempr*boltzmann # eV
        beta    = 1.0/kbT

        kappa   = 0.0 # parameter in the broadening function

        if conv
            print_conv = true
            hcg = hcg_conv(mumn, A, B, NE, fermi_Ha, freq, beta, kappa, gph, print_conv)
        else
            hcg = hcg_auto(mumn, A, B, NE, fermi_Ha, freq, beta, kappa, gph)
        end

        writedlm(name, hcg)

    end


    if dis_b
        println("Computing the population distribution.")
        dist_desc = dist_builder[1] 

        # number of total Chebyshev polynomials
        # description, temperature, lifetime, number of chebyshev polynomials,
        # number of random vectors, number of blocks, convergence, name
        desc, tempr, Gamma, N, NR, NB, NE, conv, name  = dist_builder 

        # number of Chebyshev polynomials per block
        NTx = NB
        NTy = NB

        # number of blocks
        NNL = N÷NTx + 1
        NNR = N÷NTy + 1

        # Build the Chebyshev matrix
        mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, NR)

        # Get inverse temperature
        kbT     = tempr*boltzmann # eV
        beta    = 1.0/kbT

        # Number of integration points in the energy discretization
        Nx = N*2 + 1
        Ny = Nx  
        perc1 = 100
        perc2 = 100
        Ei,Ej,opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc1, perc2)
        NEp = Nx
        Gammaf = 1
        write=false

        dist = sum_dist(opt, Ei, Ej, NE, NEp, fermi_Ha, freq, beta, Gammaf, Gamma, write)

        writedlm(name, dist)

    end

    if dos_b
        println("Computing the density of states.")

        # Description, number of energies, minimum energy in Ha, maximum energy in Ha
        # number of polynomials, number of random vectors, name of output file
        description, NE, minE_Ha, maxE_Ha, N, Nk, outname = dos_builder

        flag = 2     # flag=2 means random vectors
        perc = 100

        # Build the Chebyshev vector
        mu = doscompute_mu!(H, N, Nk, flag)

        # Resum the Chebyshev vector into the DoS
        dos = get_dos(mu, perc, minE_Ha, maxE_Ha, NE, A, B)

        # Save to file
        open(outname, "w") do io
            writedlm(io, dos)
        end
    end



    if ban_b
        println("Computing the band structure.")

        # Description and name of output file
        description, outname = band_builder

        # Calculate the bandstructure along a predefined path
        bands = get_bands(onsite, first_neighbour, second_neighbour)

        # Save to file
        open(outname, "w") do io
            writedlm(io, bands)
        end
    end
    println("Program completed successfully.")
end
