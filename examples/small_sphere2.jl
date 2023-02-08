# Configuration script

include("../src/nanoparticle.jl")

# Define the Hamiltonian
hamiltonian_builder = configure_ham_predef(shape="cube", material="gold", size="3nm")
hamiltonian_builder = configure_ham_custom(location="hamiltonian.dat")

# Select the potential to be used to generate e-h pairs
potential_builder = configure_potential_custom(location = "", frequency = "2.3 eV")
potential_builder = configure_potential_linear(frequency = "2.3eV") # use the built-in dielectric function
potential_builder = configure_potential_linear_diel(dielectric_relative = 3.1 + 1.0j, frequency = "2 ev")

# Select the output quantity: Hot carrier generation rate, or state distribution
calculation_builder = configure_hcg(Ncheb_polys = 512, Nrandom = 5)      # calculate just the hot carrier generation
calculation_builder = configure_dist(Ncheb_polys = 512, Nrandom = 5)     # calculate just the distribution
calculation_builder = configure_hcg_dist(Ncheb_polys = 512, Nrandom = 5) # calculate both the hot carrier generation and distribution


# uses the information about the Hamiltonian to print out the atomic positions 
nanoparticle_emptypot(hamiltonian_builder)

# run the whole code
nanoparticle_run(hamiltonian_builder, potential_builder, calculation_builder)

