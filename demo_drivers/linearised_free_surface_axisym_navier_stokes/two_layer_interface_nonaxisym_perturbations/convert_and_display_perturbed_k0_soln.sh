#!/usr/bin/sh

cd $@
ls perturbed_k0_soln*.dat | while read file
do
    ../../../../bin/oomph-convert.py $file
done

../../../../bin/makePvd perturbed_k0_soln perturbed_k0_soln.pvd

paraview perturbed_k0_soln.pvd

