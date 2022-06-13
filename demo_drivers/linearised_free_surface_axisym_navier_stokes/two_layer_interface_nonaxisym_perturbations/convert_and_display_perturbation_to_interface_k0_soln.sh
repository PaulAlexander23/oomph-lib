#!/usr/bin/sh

cd $@
ls perturbation_to_interface_k0_soln*.dat | while read file
do
    ../../../../bin/oomph-convert.py -p3 $file
done

../../../../bin/makePvd perturbation_to_interface_k0_soln perturbation_to_interface_k0_soln.pvd

paraview perturbation_to_interface_k0_soln.pvd

