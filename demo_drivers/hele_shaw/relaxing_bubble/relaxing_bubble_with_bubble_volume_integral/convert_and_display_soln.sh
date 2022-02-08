#!/usr/bin/sh

cd $@
ls * | while read file
do
    ../../../../../bin/oomph-convert.py -p2 $file
done

../../../../bin/makePvd soln soln.pvd

paraview soln.pvd

