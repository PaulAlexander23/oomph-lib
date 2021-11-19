#!/usr/bin/sh

~/Repositories/oomph-lib/bin/oomph-convert.py soln0.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln1.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln2.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln3.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln4.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln5.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln6.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln7.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln8.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln9.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln10.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln11.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln12.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln13.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln14.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln15.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln16.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln17.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln18.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln19.dat
~/Repositories/oomph-lib/bin/oomph-convert.py soln20.dat

~/Repositories/oomph-lib/bin/makePvd soln soln.pvd

paraview soln.pvd
