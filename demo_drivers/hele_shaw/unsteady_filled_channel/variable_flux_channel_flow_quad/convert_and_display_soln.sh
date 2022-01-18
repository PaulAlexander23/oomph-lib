#!/usr/bin/sh

cd $@
ls soln*.dat | while read file
do
   ../../../../../bin/oomph-convert.py $file
done

../../../../../bin/makePvd soln soln.pvd

paraview soln.pvd


