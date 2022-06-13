#!/usr/bin/sh

cd $@
ls perturbation_to_interface_k2_soln*.dat | while read file
do
    ../../../../bin/oomph-convert.py $file
done

../../../../bin/makePvd perturbation_to_interface_k2_soln perturbation_to_interface_k2_soln.pvd

paraview perturbation_to_interface_k2_soln.pvd

