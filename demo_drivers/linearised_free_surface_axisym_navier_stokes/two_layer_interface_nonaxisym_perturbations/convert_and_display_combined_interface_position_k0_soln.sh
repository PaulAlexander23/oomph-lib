#!/usr/bin/sh

cd $@
ls combined_interface_position_k0_soln*.dat | while read file
do
    ../../../../bin/oomph-convert.py -p2 $file
done

../../../../bin/makePvd combined_interface_position_k0_soln combined_interface_position_k0_soln.pvd

paraview combined_interface_position_k0_soln.pvd

