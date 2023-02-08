#!/bin/bash

# This script uses the constant optical matrix 'constant_phi.h5' and integrates it
# using the gw weight function along the Ep direction. This should provide a 
# function of E which is zero, then plateaus to a negative value at Ef - hw, 
# quickly plateaus to the symmetrical positive value at Ef and then goes to zero
# at Ef + hw
#
#

h5name="constant_phi.h5"
outname="dist.dat"
N1=951
N2=931
fermi=6.0
hw=2.0
kbT=0.01
debug=2

julia ../src/dist_if.jl $h5name $outname $N1 $N2 $fermi $hw $kbT $debug


