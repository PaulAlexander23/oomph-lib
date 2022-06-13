#!/usr/bin/sh

cd $@
ls base_soln*.dat | while read file
do
    ../../../../bin/oomph-convert.py $file
done

../../../../bin/makePvd base_soln base_soln.pvd

paraview base_soln.pvd

