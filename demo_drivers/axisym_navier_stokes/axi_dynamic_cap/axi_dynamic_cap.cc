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

// Solve for the steady state of a axisymmetric fluid surface in a capillary,
// under the force of gravity and a moving wall.
// Reads from a parameter file, see the default parameter file, and outputs
// bulk and surface data and restart files.

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "projectable_axisymmetric_Ttaylor_hood_elements.h"
#include "singular_axisym_dynamic_cap_problem.h"
#include "parameters.h"
#include "utility_functions.h"

using namespace std;
using namespace oomph;

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

  // Debug the jacobian
  CommandLineArgs::specify_command_line_flag("--debug-jacobian");
  CommandLineArgs::specify_command_line_flag("--pin-solid");

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
  SingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>
    problem(&parameters);
  problem.setup();
  save_dofs_types<SingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>*>(&problem, "dofs_types.dat");
  DoubleVector residuals;
  CRDoubleMatrix jacobian;
  problem.get_jacobian(residuals, jacobian);
  jacobian.sparse_indexed_output("j.dat");

  if (CommandLineArgs::command_line_flag_has_been_set("--debug-jacobian"))
  {
    debug_jacobian<SingularAxisymDynamicCapProblem<
      SingularAxisymNavierStokesElement<
        ProjectableAxisymmetricTTaylorHoodPVDElement>,
      BDF<2>>*>(&problem);
  }

  // Load in restart file
  if (parameters.restart_filename != "")
  {
    try
    {
      ifstream restart_filestream;
      restart_filestream.open(parameters.restart_filename);
      bool is_unsteady_restart = false;
      problem.read(restart_filestream, is_unsteady_restart);
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
  problem.open_trace_files(true);

  // Save a copy of the parameters
  save_parameters_to_file(parameters,
                          parameters.output_directory + "/parameters.dat");

  // Document initial condition
  problem.create_restart_file();
  problem.doc_solution();
  problem.use_fd_jacobian_for_the_bulk_augmented();

  // If the final time is zero (or less) then we are doing a steady solve,
  if (parameters.final_time <= 0)
  {
    // Solve for the steady state adapting if needed by the Z2 error estimator
    problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);

    // Create the restart file - needed before the doc solution
    problem.create_restart_file();

    // Document the solution
    problem.doc_solution();
  }
  // ...otherwise, we are doing an unsteady run
  else
  {
    // Timestep until the desired final time
    problem.timestep(parameters.time_step, parameters.final_time);
  }

  // Close the trace files
  problem.close_trace_files();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  // Return 0 to tell everyone that the program finished successfully
  return 0;
}
