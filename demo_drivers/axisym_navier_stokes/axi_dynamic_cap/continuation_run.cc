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

void normal_continuation_run(Params& parameters,
                             double& starting_step,
                             double*& continuation_param_pt)
{
  // Construct the problem
  bool has_restart = false;
  if (parameters.restart_filename != "")
  {
    std::cout << "restarting" << std::endl;
    has_restart = true;
  }
  SingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>
    problem(parameters);

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
      throw std::invalid_argument(
        "Restart filename can't be set, or opened, or read.");
    }
  }
  problem.setup();

  // Setup trace file
  problem.open_trace_files(true);

  // Save a copy of the parameters
  save_parameters_to_file(parameters,
                          parameters.output_directory + "/parameters.dat");

  // Document initial condition
  problem.create_restart_file();
  problem.doc_solution();

  double step = starting_step;
  *continuation_param_pt -= step;
  const unsigned max_iterations =
    floor(abs(parameters.final_time / parameters.time_step));
  unsigned iterations = 0;
  while (iterations < max_iterations)
  {
    // Output info
    std::cout << "Iteration: " << iterations
              << ", param: " << *continuation_param_pt << ", step: " << step
              << std::endl;

    // Update iteration number
    iterations++;

    // Continue parameter
    *continuation_param_pt += step;

    // Store current state of problem
    DoubleVector dofs;
    problem.get_dofs(dofs);
    try
    {
      // Solve for the steady state adapting if needed by the Z2 error
      // estimator
      problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);

      // Create the restart file - needed before the doc solution
      problem.create_restart_file();

      // Document the solution
      problem.doc_solution();
    }
    // If error in solving for the steady state
    catch (OomphLibException& e)
    {
      // Restore state of the problem
      problem.set_dofs(dofs);

      // Undo step
      *continuation_param_pt -= step;

      // Reduce step size
      step /= 2.0;
    }
  }

  // Close the trace files
  problem.close_trace_files();
}

void arc_continuation_run(Params& parameters,
                          double& starting_step,
                          double*& continuation_param_pt)
{
  //  // Construct the problem
  //  bool has_restart = false;
  //  if (parameters.restart_filename != "")
  //  {
  //    std::cout << "restarting" << std::endl;
  //    has_restart = true;
  //  }
  //
  //  SingularAxisymDynamicCapProblem<
  //    SingularAxisymNavierStokesElement<
  //      ProjectableAxisymmetricTTaylorHoodPVDElement>,
  //    BDF<2>>
  //    problem(Global_Physical_Params::Equilibrium_contact_angle,
  //    has_restart);
  //
  //  // Load in restart file
  //  if (parameters.restart_filename != "")
  //  {
  //    try
  //    {
  //      ifstream restart_filestream;
  //      restart_filestream.open(parameters.restart_filename);
  //      bool is_unsteady_restart = false;
  //      problem.read(restart_filestream, is_unsteady_restart);
  //      restart_filestream.close();
  //    }
  //    catch (exception& e)
  //    {
  //      throw std::invalid_argument(
  //        "Restart filename can't be set, or opened, or read.");
  //    }
  //  }
  //
  //  problem.set_contact_angle(
  //    Global_Physical_Params::Equilibrium_contact_angle);
  //  problem.set_bond_number(Global_Physical_Params::Bo);
  //  problem.set_capillary_number(Global_Physical_Params::Ca);
  //  problem.set_reynolds_number(Global_Physical_Params::Re);
  //
  //  // Set maximum number of mesh adaptations per solve
  //  problem.set_max_adapt(parameters.max_adapt);
  //
  //  // Set output directory
  //  problem.set_directory(parameters.dir_name);
  //
  //  // Setup trace file
  //  problem.open_trace_files(true);
  //
  //  ofstream parameters_filestream(
  //    (parameters.dir_name + "/parameters.dat").c_str());
  //  parameters.doc(parameters_filestream);
  //  parameters_filestream.close();
  //
  //  // Document initial condition
  //  problem.create_restart_file();
  //  problem.doc_solution();
  //
  //  // Solve for the steady state adapting if needed by the Z2 error estimator
  //  problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  //
  //  // Create the restart file - needed before the doc solution
  //  problem.create_restart_file();
  //
  //  // Document the solution
  //  problem.doc_solution();
  //
  //  // Set the tracking parameter
  //  problem.set_analytic_dparameter(
  //    problem.reynolds_number_inverse_froude_number_pt());
  //  // problem.set_analytic_dparameter(&Slip_Params::wall_velocity);
  //
  //  double* parameter_pt = 0;
  //  if (continuation_param_pt == &Global_Physical_Params::Bo)
  //  {
  //    parameter_pt = problem.reynolds_number_inverse_froude_number_pt();
  //  }
  //  else if (continuation_param_pt == &Slip_Params::wall_velocity)
  //  {
  //    parameter_pt = &Slip_Params::wall_velocity;
  //  }
  //  else if (continuation_param_pt ==
  //           &Global_Physical_Params::Equilibrium_contact_angle)
  //  {
  //    parameter_pt = problem.get_contact_angle_pt();
  //  }
  //  else
  //  {
  //    throw std::runtime_error("Not implemented yet.");
  //  }
  //
  //  double ds = starting_step;
  //  const unsigned number_of_steps = floor(abs(parameters.ft /
  //  parameters.dt)); for (unsigned n = 1; n < number_of_steps; n++)
  //  {
  //    problem.arc_length_step_solve(parameter_pt, ds);
  //    problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  //
  //    problem.create_restart_file();
  //    problem.doc_solution();
  //
  //    // Adapt and solve the problem by the number of intervals between adapts
  //    // parameter.
  //    // if (n % Mesh_Control_Params::interval_between_adapts ==
  //    //    Mesh_Control_Params::interval_between_adapts - 1)
  //    //{
  //    //  // Solve for the steady state adapting if needed by the Z2 error
  //    //  estimator problem.reset_arc_length_parameters();
  //    //  problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  //    //}
  //  }
  //
  //  // Close the trace files
  //  problem.close_trace_files();
}

void bond_height_control_continuation_run(Params& parameters,
                                          double& starting_step,
                                          double*& continuation_param_pt);

void ca_height_control_continuation_run(Params& parameters,
                                        double& starting_step,
                                        double*& continuation_param_pt);

void height_control_continuation_run(Params& parameters,
                                     double& starting_step,
                                     double*& continuation_param_pt)
{
  //  if (continuation_param_pt == &Global_Physical_Params::Bo)
  //  {
  //    bond_height_control_continuation_run(
  //      parameters, starting_step, continuation_param_pt);
  //  }
  //  else if (continuation_param_pt == &Slip_Params::wall_velocity)
  //  {
  //    ca_height_control_continuation_run(
  //      parameters, starting_step, continuation_param_pt);
  //  }
}

void bond_height_control_continuation_run(Params& parameters,
                                          double& starting_step,
                                          double*& continuation_param_pt)
{
  //  // Construct the problem
  //  bool has_restart = false;
  //  if (parameters.restart_filename != "")
  //  {
  //    std::cout << "restarting" << std::endl;
  //    has_restart = true;
  //  }
  //  BoHeightControlSingularAxisymDynamicCapProblem<
  //    SingularAxisymNavierStokesElement<
  //      ProjectableAxisymmetricTTaylorHoodPVDElement>,
  //    BDF<2>>
  //    problem(Global_Physical_Params::Equilibrium_contact_angle,
  //    has_restart);
  //
  //  // Load in restart file
  //  if (parameters.restart_filename != "")
  //  {
  //    try
  //    {
  //      ifstream restart_filestream;
  //      restart_filestream.open(parameters.restart_filename);
  //      bool is_unsteady_restart = false;
  //      problem.read(restart_filestream, is_unsteady_restart);
  //      restart_filestream.close();
  //    }
  //    catch (exception& e)
  //    {
  //      throw std::invalid_argument(
  //        "Restart filename can't be set, or opened, or read.");
  //    }
  //  }
  //
  //  problem.set_contact_angle(
  //    Global_Physical_Params::Equilibrium_contact_angle);
  //  problem.set_bond_number(Global_Physical_Params::Bo);
  //  problem.set_capillary_number(Global_Physical_Params::Ca);
  //  problem.set_reynolds_number(Global_Physical_Params::Re);
  //  problem.set_max_adapt(parameters.max_adapt);
  //  problem.set_directory(parameters.dir_name);
  //  problem.open_trace_files(true);
  //
  //  ofstream parameters_filestream(
  //    (parameters.dir_name + "/parameters.dat").c_str());
  //  parameters.doc(parameters_filestream);
  //  parameters_filestream.close();
  //
  //  // Document initial condition
  //  problem.create_restart_file();
  //  problem.doc_solution();
  //
  //  // Solve for the steady state adapting if needed by the Z2 error estimator
  //  problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  //
  //  // Document the solution
  //  problem.create_restart_file();
  //  problem.doc_solution();
  //
  //  if (continuation_param_pt == &Global_Physical_Params::Bo)
  //  {
  //    problem.set_height_step_parameter_to_bond_number();
  //  }
  //  else if (continuation_param_pt == &Slip_Params::wall_velocity)
  //  {
  //    problem.set_height_step_parameter_to_wall_velocity();
  //  }
  //  else
  //  {
  //    throw std::runtime_error("Not implemented yet.");
  //  }
  //
  //  double ds = starting_step;
  //  const unsigned number_of_steps = floor(abs(parameters.ft /
  //  parameters.dt)); for (unsigned n = 1; n < number_of_steps; n++)
  //  {
  //    TerminateHelper::setup();
  //    problem.height_step_solve(ds);
  //    problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  //
  //    problem.create_restart_file();
  //    problem.doc_solution();
  //  }
  //
  //  // Close the trace files
  //  problem.close_trace_files();
  //  TerminateHelper::clean_up_memory();
}

void ca_height_control_continuation_run(Params& parameters,
                                        double& starting_step,
                                        double*& continuation_param_pt)
{
  //  // Construct the problem
  //  bool has_restart = false;
  //  if (parameters.restart_filename != "")
  //  {
  //    std::cout << "restarting" << std::endl;
  //    has_restart = true;
  //  }
  //  CaHeightControlSingularAxisymDynamicCapProblem<
  //    SingularAxisymNavierStokesElement<
  //      ProjectableAxisymmetricTTaylorHoodPVDElement>,
  //    BDF<2>>
  //    problem(Global_Physical_Params::Equilibrium_contact_angle,
  //    has_restart);
  //
  //  // Load in restart file
  //  if (parameters.restart_filename != "")
  //  {
  //    try
  //    {
  //      ifstream restart_filestream;
  //      restart_filestream.open(parameters.restart_filename);
  //      bool is_unsteady_restart = false;
  //      problem.read(restart_filestream, is_unsteady_restart);
  //      restart_filestream.close();
  //    }
  //    catch (exception& e)
  //    {
  //      throw std::invalid_argument(
  //        "Restart filename can't be set, or opened, or read.");
  //    }
  //  }
  //
  //  problem.set_contact_angle(
  //    Global_Physical_Params::Equilibrium_contact_angle);
  //  problem.set_bond_number(Global_Physical_Params::Bo);
  //  problem.set_capillary_number(Global_Physical_Params::Ca);
  //  problem.set_reynolds_number(Global_Physical_Params::Re);
  //  problem.set_max_adapt(parameters.max_adapt);
  //  problem.set_directory(parameters.dir_name);
  //  problem.open_trace_files(true);
  //
  //  ofstream parameters_filestream(
  //    (parameters.dir_name + "/parameters.dat").c_str());
  //  parameters.doc(parameters_filestream);
  //  parameters_filestream.close();
  //
  //  // Document initial condition
  //  problem.create_restart_file();
  //  problem.doc_solution();
  //
  //  // Solve for the steady state adapting if needed by the Z2 error estimator
  //  problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  //
  //  // Document the solution
  //  problem.create_restart_file();
  //  problem.doc_solution();
  //
  //  if (continuation_param_pt == &Global_Physical_Params::Bo)
  //  {
  //    problem.set_height_step_parameter_to_bond_number();
  //  }
  //  else if (continuation_param_pt == &Slip_Params::wall_velocity)
  //  {
  //    problem.set_height_step_parameter_to_wall_velocity();
  //  }
  //  else
  //  {
  //    throw std::runtime_error("Not implemented yet.");
  //  }
  //
  //  double ds = starting_step;
  //  const unsigned number_of_steps = floor(abs(parameters.ft /
  //  parameters.dt)); for (unsigned n = 1; n < number_of_steps; n++)
  //  {
  //    TerminateHelper::setup();
  //    try
  //    {
  //      problem.height_step_solve(ds);
  //      problem.steady_newton_solve_adapt_if_needed(parameters.max_adapt);
  //    }
  //    catch (exception& e)
  //    {
  //      ds /= 2;
  //    }
  //
  //    problem.create_restart_file();
  //    problem.doc_solution();
  //  }
  //
  //  // Close the trace files
  //  problem.close_trace_files();
  //  TerminateHelper::clean_up_memory();
}

int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  // Setup mpi but don't make a copy of mpi_comm_world because
  // mumps wants to work with the real thing.
  bool make_copy_of_mpi_comm_world = false;
  MPI_Helpers::init(argc, argv, make_copy_of_mpi_comm_world);
#endif

  // Parse the arguments
  Arguments args = parse_arguments(argc, argv);

  // Problem parameters
  Params parameters = create_parameters_from_file(args.parameters_filename);

  double* continuation_param_pt = 0;
  // Continuation parameter pointer
  switch (args.continuation_param)
  {
    case ContinuationParameter::Bo:
      continuation_param_pt = &parameters.reynolds_inverse_froude_number;
      break;
    case ContinuationParameter::Ca:
      continuation_param_pt = &parameters.wall_velocity;
      break;
    case ContinuationParameter::Angle:
      continuation_param_pt = &parameters.contact_angle;
      break;
    case ContinuationParameter::SlipLength:
      continuation_param_pt = &parameters.slip_length;
      break;
  }


  if (args.has_arc_continuation)
  {
    arc_continuation_run(parameters, args.starting_step, continuation_param_pt);
  }
  else if (args.has_height_control_continuation)
  {
    height_control_continuation_run(
      parameters, args.starting_step, continuation_param_pt);
  }
  else
  {
    normal_continuation_run(
      parameters, args.starting_step, continuation_param_pt);
  }

// Finalise MPI after all computations are complete
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  return 0;
}
