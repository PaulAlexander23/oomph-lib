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

// OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "singular_axisym_navier_stokes_elements.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"

// Local include files
#include "hijacked_projectable_axisymmteric_Ttaylor_hood_elements.h"
#include "bo_height_control_singular_axisym_dynamic_cap_problem.h"
#include "parameters.h"

using namespace std;
using namespace oomph;

int main(int argc, char** argv)
{
  // Check number of arguments
  int number_of_arguments = argc - 1;
  if (number_of_arguments != 3)
  {
    cout << "Wrong number of arguments." << std::endl;
    cout << "Please provide either --wall_velocity, --Bo or --angle with a "
            "starting step ds."
         << std::endl;
    return 1;
  }

#ifdef OOMPH_HAS_MPI
  // Setup mpi but don't make a copy of mpi_comm_world because
  // mumps wants to work with the real thing.
  bool make_copy_of_mpi_comm_world = false;
  MPI_Helpers::init(argc, argv, make_copy_of_mpi_comm_world);
#endif

  // Problem parameters
  Parameters parameters;
  read_parameters_from_file(argv[3], parameters);

  // Construct the problem
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    cout << "restarting" << endl;
    has_restart = true;
  }
  BoHeightControlSingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      HijackedProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>
    problem(Global_Physical_Parameters::Equilibrium_contact_angle, has_restart);

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
      cout << "Restart filename can't be set, or opened, or read." << endl;
      cout << "File: " << parameters.restart_filename << endl;
      return 1;
    }
  }

  problem.set_contact_angle(
    Global_Physical_Parameters::Equilibrium_contact_angle);
  problem.set_bond_number(Global_Physical_Parameters::Bo);
  problem.set_capillary_number(Global_Physical_Parameters::Ca);
  problem.set_reynolds_number(Global_Physical_Parameters::Re);
  problem.set_max_adapt(parameters.max_adapt);
  problem.set_directory(parameters.dir_name);
  problem.open_trace_files(true);

  ofstream parameters_filestream(
    (parameters.dir_name + "/parameters.dat").c_str());
  parameters.doc(parameters_filestream);
  parameters_filestream.close();

  // Document initial condition
  problem.create_restart_file();
  problem.doc_solution();

  // Solve for the steady state adapting if needed by the Z2 error estimator
  problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);

  // Document the solution
  problem.create_restart_file();
  problem.doc_solution();

  std::ofstream output_stream;
  output_stream.open("dofs.txt");
  problem.describe_dofs(output_stream);
  output_stream.close();

  try
  {
    if (string(argv[1]) == "--Bo")
    {
      problem.set_height_step_parameter_to_bond_number();
    }
    else if (string(argv[1]) == "--wall_velocity")
    {
      problem.set_height_step_parameter_to_wall_velocity();
    }
    else
    {
      cout << "Not implemented yet for this input. Arg in: " << argv[1] << endl;
    }
  }
  catch (exception& e)
  {
    cout << "Couldn't set the continuation parameter. Arg in: " << argv[1]
         << endl;
  }

  double ds = 0;
  try
  {
    ds = atof(argv[2]);
  }
  catch (exception& e)
  {
    cout << "Couldn't set arc length step. Arg in: " << argv[2] << endl;
  }

  const unsigned number_of_steps = 100;
  for (unsigned n = 0; n < number_of_steps; n++)
  {
    TerminateHelper::setup();
    problem.height_step_solve(ds);

    problem.create_restart_file();
    problem.doc_solution();

    // Adapt and solve the problem by the number of intervals between adapts
    // parameter.
    if (n % Mesh_Control_Parameters::interval_between_adapts ==
        Mesh_Control_Parameters::interval_between_adapts - 1)
    {
      // Solve for the steady state adapting if needed by the Z2 error estimator
      problem.reset_arc_length_parameters();
      if (string(argv[1]) == "--Bo")
      {
        problem.adapt();
      }
      problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
    }

    // If the contact angle error gets too large
    if (problem.get_contact_angle_error() > 1.0)
    {
      // Stop the loop
      break;
    }
  }

  // Close the trace files
  problem.close_trace_files();
  TerminateHelper::clean_up_memory();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  return 0;
}