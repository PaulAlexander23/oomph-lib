# Hele-Shaw bubble

Extracted from the suplementary material for the article "The influence of invariant solutions on the transient behaviour of an air bubble in a Hele-Shaw channel." published in Proc. R. Soc. A, 2019 by J. S. Keeler et al.

## Installation

1. Clone or download oomph-lib from the Git repository: https://github.com/oomph-lib/oomph-lib
2. Download trillinos. See http://oomph-lib.maths.man.ac.uk/doc/html/index.html
3. Place trillinos tar file in `oomph-lib/external_distributions/trillinos/`
6. Run `./autogen.sh` in parent directory and follow prompts
7. Update `OOMPH-LIB_INSTALL_LOCATION` in Makefile of this repository
8. Run `make`

## Testing

Testing is easily run via,
```
make test
```

## Usage

* `./bubble_steady` - Calculates a steady solution and continues a solution branch from a given restart file given as a command line argument.
* `./bubble_unsteady` - Starts from an arbitrary initial condition and time-steps. Choose the parameters for a given initial condition when prompted.
* `./bubble_steady_weakly_nonlinear` - Calculates the Landau coefficients from the restart_H1 and restart_H2 files. For converged results perform a bisection but this will take a long time...

Note: These all require command line inputs. Add a `--help` to display the options.

### Restart files

The restart files, which can be found at https://figshare.com/s/254feb5ce70a340b0d57?file=14702084, can be used to read in a solution for bubble_steady and bubble_steady_weakly_nonlinear.
e.g. on the command line,
```
./bubble_steady -n 4 -f restart_S1.dat
```
will restart a solution branch from the S1 solution (see reference_figure.pdf in the suplementary material for a correspondence to each solution) whilst
```
./bubble_steady_weakly_nonlinear -d 1 -f restart_H1.dat
```
will calculate the Landau coefficients for the H1 point.
