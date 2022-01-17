#!/usr/bin/sh

cd RESLT
../../../../bin/oomph-convert.py soln.dat
paraview soln.vtu
