

function nanoparticle_run(ham_builder, pot_builder, cal_builder)


    # Obtain the description for each section of the workflow
    pot_desc = pot_builder[1]
    ham_desc = ham_builder[1]
    cal_desc = cal_builder[1]




    # Get Hamiltonian from the predefined set
    if ham_desc == "predefined"
        shape = ham_desc[2]
        mater = ham_desc[3] 
        size  = ham_desc[4]

        H, R



    # Get Hamiltonian sparse matrix from file
    elseif ham_desc == "custom"
        ham_fname = ham_builder[2]
        # H = ()

    else 
        println("Hamiltonian description not supported.")
    end


    # Dielectric sphere, use the dielectric of the material
    if pot_desc == "linear"



    # Dielectric sphere, use custom dielectric constant
    elseif pot_desc == "linear_diel"
        
    # Reading from file
    elseif pot_desc == "custom"
        pot_fname = pot_builder[2]

    else
        print("Potential description not supported.")
    end



end
