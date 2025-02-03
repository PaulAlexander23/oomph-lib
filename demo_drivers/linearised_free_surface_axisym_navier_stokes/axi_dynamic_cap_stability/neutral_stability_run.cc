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

#include "algorithm"

// OOMPH-LIB include files
#include "generic.h"
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

enum ContinuationParameter
{
  Bo,
  Ca,
  Angle,
  SlipLength
};

struct Arguments
{
  std::string parameters_filename = "";
  double starting_step = 0.0;
  ContinuationParameter continuation_param = ContinuationParameter::Bo;
  bool has_arc_continuation = false;
  bool has_height_control_continuation = false;
};

Arguments parse_arguments(const int& argc, char** const& argv)
{
  Arguments args;

  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
    "--parameters", &args.parameters_filename, "Required: parameter filename");

  double bo_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--Bo",
    &bo_initial_step,
    "Optional: Continue in Bond number with initial step");

  double ca_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--wall_velocity",
    &ca_initial_step,
    "Optional: Continue in Capillary number with initial step");

  double angle_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--angle",
    &angle_initial_step,
    "Optional: Continue in angle with initial step");

  double slip_length_initial_step = 0.0;
  CommandLineArgs::specify_command_line_flag(
    "--slip_length",
    &slip_length_initial_step,
    "Optional: Continue in slip length with initial step");

  CommandLineArgs::specify_command_line_flag("--arc",
                                             "Optional: Use arc continuation");

  CommandLineArgs::specify_command_line_flag(
    "--height_control", "Optional: Use height control continuation");


  // Parse and assign command line arguments
  bool has_unrecognised_arg = false;
  CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);
  if (has_unrecognised_arg)
  {
    throw std::invalid_argument("Unrecognised args.");
  }
  if (args.parameters_filename == "")
  {
    throw std::invalid_argument("Parameter file not set.");
  }

  // Set continuation variables
  try
  {
    if (bo_initial_step != 0)
    {
      args.starting_step = bo_initial_step;
      args.continuation_param = ContinuationParameter::Bo;
    }
    else if (angle_initial_step != 0)
    {
      args.starting_step = angle_initial_step;
      args.continuation_param = ContinuationParameter::Angle;
    }
    else if (ca_initial_step != 0)
    {
      args.starting_step = ca_initial_step;
      args.continuation_param = ContinuationParameter::Ca;
    }
    else if (slip_length_initial_step != 0)
    {
      args.starting_step = slip_length_initial_step;
      args.continuation_param = ContinuationParameter::SlipLength;
    }
    else
    {
      throw std::invalid_argument("No continuation parameter specified.");
    }
  }
  catch (exception& e)
  {
    throw std::invalid_argument("Continuation parameter can't be set.");
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--arc"))
  {
    args.has_arc_continuation = true;
  }
  if (CommandLineArgs::command_line_flag_has_been_set("--height_control"))
  {
    args.has_height_control_continuation = true;
  }
  return args;
}

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

  Arguments args = parse_arguments(argc, argv);

  // Problem parameters
  Params parameters = create_parameters_from_file(args.parameters_filename);

  double* continuation_param_pt = 0;
  // Continuation parameter pointer
  switch (args.continuation_param)
  {
    case ContinuationParameter::Bo:
      continuation_param_pt = parameters.reynolds_inverse_froude_number_pt;
      break;
    case ContinuationParameter::Ca:
      continuation_param_pt = parameters.wall_velocity_pt;
      break;
    case ContinuationParameter::Angle:
      continuation_param_pt = &parameters.contact_angle;
      break;
    case ContinuationParameter::SlipLength:
      continuation_param_pt = &parameters.slip_length;
      break;
  }

  // Construct the base problem
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    std::cout << "restarting" << std::endl;
    has_restart = true;
  }
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
    catch (OomphLibError& e)
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

  //====================================================================

  // Solve steady problem
  base_problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Output result
  base_problem.create_restart_file();
  base_problem.doc_solution();

  //====================================================================
  // Solve the eigenvalue problem
  Vector<std::complex<double>> eigenvalues(1, 0.0);


  // Create the linear stability problem
  typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                     PERTURBED_ELEMENT,
                                     TIMESTEPPER>
    initial_perturbed_problem(base_problem.bulk_mesh_pt(),
                              base_problem.free_surface_mesh_pt(),
                              base_problem.slip_surface_mesh_pt(),
                              &parameters);

  initial_perturbed_problem.assign_initial_values_impulsive();

  // Document the solution before the solve for testing
  initial_perturbed_problem.doc_solution();
  initial_perturbed_problem.steady_newton_solve();

  if (parameters.azimuthal_mode_number != 0)
  {
    initial_perturbed_problem.make_unsteady();
  }

  eigenvalues =
    initial_perturbed_problem.solve_n_most_unstable_eigensolutions(1);

  // Our new residual is the real part of the eigenvalue
  double residual = real(eigenvalues[0]);

  const double tolerance = 1e-8;
  bool has_converged = false;
  if (abs(residual) < tolerance)
  {
    has_converged = true;
  }

  // Newton iteration until we converge to the neutral stability point
  //====================================================================

  unsigned n_iterations = 0;
  const unsigned max_n_iterations =
    floor(abs(parameters.final_time / parameters.time_step));
  // const double step_tolerance = 1e-3;
  double step_param = args.starting_step;
  double old_param;
  double current_param = *continuation_param_pt;
  double old_residual;
  while (!has_converged && n_iterations < max_n_iterations)
  {
    TerminateHelper::setup();
    std::cout << "-------------" << std::endl;
    std::cout << "Start of loop" << std::endl;
    std::cout << "-------------" << std::endl;
    std::cout << std::endl;
    std::cout << n_iterations << ", " << current_param << "," << residual << ","
              << current_param + step_param << std::endl;
    old_param = current_param;
    current_param += step_param;
    *continuation_param_pt = current_param;


    // Solve steady problem
    bool has_base_state = false;
    unsigned base_state_iterations = 0;
    const unsigned max_base_state_iterations = 8;
    while (!has_base_state && base_state_iterations < max_base_state_iterations)
    {
      DoubleVector dofs;
      base_problem.get_dofs(dofs);
      try
      {
        // Solve steady problem
        int exit_flag = base_problem.steady_newton_solve_adapt_if_needed(
          parameters.max_adapt);

        // Output result
        base_problem.create_restart_file();
        base_problem.doc_solution();

        if (exit_flag >= 0)
        {
          has_base_state = true;
        }
      }
      catch (OomphLibError& e)
      {
        std::cout << "Caught exception" << std::endl;
        std::cout << "Resetting problem" << std::endl;
        base_problem.set_dofs(dofs);
        current_param -= step_param;
        step_param /= 3.0;
        current_param += step_param;
        *continuation_param_pt = current_param;
        std::cout << "Reducing step size.";
        std::cout << "Number of attempts: " << base_state_iterations;
        std::cout << ", Step size: " << step_param;
        std::cout << ", Target wall velocity: " << current_param << std::endl;
      }
      base_state_iterations++;
    }

    if (!has_base_state)
    {
      std::cout << "WARNING: Base state not found." << std::endl;
      break;
    }

    //====================================================================

    // Create the linear stability problem
    PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                       PERTURBED_ELEMENT,
                                       TIMESTEPPER>
      perturbed_problem(base_problem.bulk_mesh_pt(),
                        base_problem.free_surface_mesh_pt(),
                        base_problem.slip_surface_mesh_pt(),
                        &parameters);

    // Document the solution before the solve for testing
    perturbed_problem.doc_solution();

    perturbed_problem.steady_newton_solve();
    if (parameters.azimuthal_mode_number != 0)
    {
      perturbed_problem.make_unsteady();
    }

    eigenvalues = perturbed_problem.solve_n_most_unstable_eigensolutions(1);

    old_residual = residual;
    residual = real(eigenvalues[0]);

    if (abs(residual) < tolerance)
    {
      has_converged = true;
      perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
    }
    else
    {
      // Newton iteration, using fd to approximate the derivative
      double deriv = (residual - old_residual) / (current_param - old_param);
      const double relaxation = 1.0;
      step_param = -relaxation * residual / (deriv);
      if (n_iterations == (max_n_iterations - 1))
      {
        perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
      }
    }

    n_iterations++;
  }

  if (!has_converged)
  {
    std::cout << "WARNING: Critical Bond number not converged." << std::endl;
  }

  // Close the trace files
  base_problem.close_trace_files();

  TerminateHelper::clean_up_memory();

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main
