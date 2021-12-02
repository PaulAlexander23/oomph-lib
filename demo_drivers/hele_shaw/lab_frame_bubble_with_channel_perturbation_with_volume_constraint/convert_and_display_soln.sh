#!/usr/bin/sh

cd $1
ls soln_*.dat | while read file
do
    ../../../../bin/oomph-convert.py $file
done

../../../../bin/makePvd soln soln.pvd

paraview soln.pvd

