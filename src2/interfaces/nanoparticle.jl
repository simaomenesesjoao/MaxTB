
boltzmann = 8.617333e-5      # [eV/Kelvin] Boltzmann constant

function nanoparticle_run(ham_builder, pot_builder, cal_builder)
    # Perform the complete workflow to obtain the hot carrier generation/distribution 


    # Obtain the description for each section of the workflow
    ham_desc = ham_builder[1]
    pot_desc = pot_builder[1]
    cal_desc = cal_builder[1]


    # Check compatibility
    if ham_desc != "predefined" and pot_desc == "linear" 
        println("Requesting a linear potential without specifying the dielectric function requires knowledge of the material, which is defined in predefined Hamiltonians. Exiting\n")
    end


    # Get Hamiltonian from the predefined set
    if ham_desc == "predefined"
        A = ham_builder[2]
        B = ham_builder[3]
        fermi = ham_builder[4]

        shape = ham_builder[5]
        mater = ham_builder[6] 
        size  = ham_builder[7]

        # Generate list of atomic positions
        Elist, Edict, R = generate_shape(size, shape)

        # Use list of atomic positions to determine the Hamiltonian
        H = slater_koster_FCC(Elist, Edict, mater)


    # Get Hamiltonian sparse matrix from file
    elseif ham_desc == "custom"
        ham_fname = ham_builder[2]
        # H, A, B, fermi = get_hami(ham_fname)

    else 
        println("Requested Hamiltonian description is not supported.")
    end



    # Frequency of the optical potential. Always required
    freq = pot_builder[2]

    # Dielectric sphere, use the dielectric of the material: pot_builder = [pot_desc, freq, eps_m]
    if pot_desc == "linear"
        eps_m = pot_builder[3] # dielectric of the environment
        eps   = get_dielectric(freq, mater)
        Phi   = potential_sphere(R, Edict, eps, eps_m)

    # Dielectric sphere, use custom dielectric constant: pot_builder = [pot_desc, freq, eps_m, eps]
    elseif pot_desc == "linear_diel"
        eps_m = pot_builder[3] # dielectric of the environment
        eps   = pot_builder[4]
        Phi   = potential_sphere(R, Edict, eps, eps_m, Norbitals = 9)

    # Reading from file: pot_builder = [pot_desc, freq, filename]
    elseif pot_desc == "custom"
        pot_fname = pot_builder[3]
        Phi = comsol_read(pot_fname)

    else
        print("Potential description not supported.")
    end




    cal_type = cal_builder[2] # HCG/dist/both
    cal_conv = cal_builder[3] # get convergence analysis

    # Simple calculation cal_builder = [cal_desc, cal_type, N]
    if cal_desc == "simple"
        # number of total Chebyshev polynomials
        N = cal_builder[4]

        # number of Chebyshev polynomials per block
        NTx = 10
        NTy = 10

        # number of blocks
        NNL = N÷NTx + 1
        NNR = N÷NTy + 1
        Nk  = 5 # number of random vectors

        # Build the Chebyshev matrix
        mumn = compute_mumn!(H, Phi, NNL, NNR, NTx, NTy, Nk)

        # Number of integration points in the energy discretization
        Nx = N*2 
        Ny = Nx  

        # Percentage of Chebyshev polynomials to keep
        perc1 = 100
        perc2 = perc1

        # Transform the Chebyshev matrix into the energy-resolved optical matrix
        Ei,Ej,opt = resum_mu_dd(mumn, Nx, Ny, A, B, perc1, perc2)

        NE  = N*3 # Number of energy points in the integration
        NEp = N*3 # Number of output energy points

        # Temperature
        tempr   = 298
        kbT     = tempr*boltzmann # eV
        beta    = 1.0/kbT

        write  = false # write debug information
        Gammaf = 1 # prefactor to the relaxation function 

        # Relaxation
        function Gamma(E,Ep)
            return 0.1
        end

        Ei *= HaeV
        Ej *= HaeV

        dist = sum_dist(opt, Ei, Ej, NE, NEp, fermi, freq, beta, Gammaf, Gamma, write)

        kappa = 0.0
        gph   = 0.06/HaeV            # [Ha] broadening
        hcg  = sum_hcg( opt, Ei, Ej, NE, NEp, fermi, freq, beta, kappa, gph, write)

end
