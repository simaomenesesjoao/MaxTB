#!/bin/bash

h5name="moments.h5"
outname="constant_phi.h5"

A=0.1
B=1.0
N1=152
N2=189
perc1=30
perc2=40

julia ../src/sumcheb_if.jl $h5name $outname $A $B $N1 $N2 $perc1 $perc2 
