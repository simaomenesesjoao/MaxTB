#!/bin/bash

# parameters to get the atomic positions
shape=sphere  # shape to be considered
rescale=false # do not rescale the results
R=2.0         # radius of the sphere [nm]

# parameters to get the potential
freq=2.0     # frequency [eV]

# parameters to get the Chebyshev matrix
NR=5        # number of random vectors
M=200        # number of Chebyshev moments

# parameters to get the hot carrier generation rate
Fermi=0.2595 # Fermi energy [Ha]
NE=1000      # number of energies to be used in the integration
perc=100     # percentage of Chebyshev polynomials to be kept

# This calculation requires around 600MB of RAM
# 1M atoms uses around 32GB of RAM with blocks of size 30

echo "script: Generating positions for the shape: $shape with dimensions $R"
julia ../src/generate_positions.jl positions.dat $shape $R $rescale

echo -e "\nscript: Determining potential. This step can also be done with Comsol"
julia ../src/potential_sphere.jl positions.dat potential.dat $freq

echo -e "\nscript: Calculating Chebyshev moments"
julia ../src/moments.jl potential.dat output.h5 $shape $R $M $NR

echo -e "\nscript: Calculating hot carrier generation rate"
julia ../src/process_moments.jl output.h5 rates.dat $freq $Fermi $NE $perc
