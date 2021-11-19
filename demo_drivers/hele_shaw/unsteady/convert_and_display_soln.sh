#!/usr/bin/sh

cd RESLT
../../../../bin/oomph-convert.py soln0.dat
../../../../bin/oomph-convert.py soln1.dat
../../../../bin/oomph-convert.py soln2.dat
../../../../bin/oomph-convert.py soln3.dat
../../../../bin/oomph-convert.py soln4.dat
../../../../bin/oomph-convert.py soln5.dat
../../../../bin/oomph-convert.py soln6.dat
../../../../bin/oomph-convert.py soln7.dat
../../../../bin/oomph-convert.py soln8.dat
../../../../bin/oomph-convert.py soln9.dat
../../../../bin/oomph-convert.py soln10.dat

../../../../bin/makePvd soln soln.pvd

paraview soln.pvd

