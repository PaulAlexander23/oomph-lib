// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Driver code for an axisymmetric free-surface hydrostatics problem.
// The system consists of a layer of fluid
// in a domain of height 1 and radius 0.5.
// The program solves for the interface position as the contact angle
// at the wall, alpha, decreases from pi/2. The resulting shapes should all be
// spherical shells and the pressure jump across the interface should be
// 2 cos(alpha)/0.5 = 4 cos(alpha)/Ca.

// STD includes
#include <algorithm>

// OOMPH-LIB include files
#include "generic.h"
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"
#include "linearised_axisym_navier_stokes.h"

// Other demo includes
#include "../../axisym_navier_stokes/axi_dynamic_cap/projectable_axisymmetric_Ttaylor_hood_elements.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/singular_axisym_dynamic_cap_problem.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/parameters.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/utility_functions.h"

// Local includes
#include "linearised_axisymmetric_fluid_interface_elements.h"
#include "decomposed_linear_elasticity_elements.h"
#include "linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_my_linear_element.h"
#include "perturbed_linear_stability_cap_problem.h"

using namespace std;
using namespace oomph;


//===start_of_main=======================================================
/// Main driver: Build problem and initiate parameter study
//======================================================================
int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  // Setup mpi but don't make a copy of mpi_comm_world because
  // mumps wants to work with the real thing.
  bool make_copy_of_mpi_comm_world = false;
  MPI_Helpers::init(argc, argv, make_copy_of_mpi_comm_world);
#endif

  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Parameter file
  std::string parameters_filename = "default_parameters.dat";
  CommandLineArgs::specify_command_line_flag("--parameters",
                                             &parameters_filename);

  // Parse command line
  const bool throw_exception_if_unrecognised_flags = true;
  CommandLineArgs::parse_and_assign(throw_exception_if_unrecognised_flags);

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // Problem parameters
  Params parameters = create_parameters_from_file(parameters_filename);

  // Construct the problem
  typedef SingularAxisymNavierStokesElement<
    ProjectableAxisymmetricTTaylorHoodPVDElement>
    BASE_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> base_problem(
    &parameters);

  // Load in restart file
  if (parameters.restart_filename != "")
  {
    try
    {
      ifstream restart_filestream;
      restart_filestream.open(parameters.restart_filename);
      bool is_unsteady_restart = false;
      base_problem.read(restart_filestream, is_unsteady_restart);
      restart_filestream.close();
    }
    catch (exception& e)
    {
      std::cout << "Restart filename can't be set, or opened, or read."
                << std::endl;
      std::cout << "File: " << parameters.restart_filename << std::endl;
      return 1;
    }
  }

  // Setup trace file
  base_problem.open_trace_files(true);

  // Save a copy of the parameters
  save_parameters_to_file(parameters,
                          parameters.output_directory + "/parameters.dat");


  // Output the initial condition
  base_problem.create_restart_file();
  base_problem.doc_solution();
  base_problem.use_fd_jacobian_for_the_bulk_augmented();

  // Solve steady problem
  base_problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Output result
  base_problem.create_restart_file();
  base_problem.doc_solution();

  // Close the trace files
  base_problem.close_trace_files();

  // Create the linear stability problem
  typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                     PERTURBED_ELEMENT,
                                     TIMESTEPPER>
    perturbed_problem(base_problem.bulk_mesh_pt(),
                      base_problem.free_surface_mesh_pt(),
                      base_problem.slip_surface_mesh_pt(),
                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();

  // Document the solution before the solve for testing
  perturbed_problem.doc_solution();

  // Solve
  // eigenproblem
  if (parameters.azimuthal_mode_number != 0)
  {
    perturbed_problem.make_unsteady();
  }
  perturbed_problem.pin_horizontal_mesh_deformation();
  perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(8);

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main
