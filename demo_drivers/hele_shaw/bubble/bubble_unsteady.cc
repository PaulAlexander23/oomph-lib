// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//           Version 0.85. June 9, 2008.
// LIC//
// LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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

// Generic routines
#include "generic.h"
#include "hele_shaw.h"
#include "meshes.h"

// The equations
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

using namespace std;
using namespace oomph;

#include "problem_parameter.h"
//#include "hele_shaw_interface_elements_with_integrals.h"
#include "custom_hele_shaw_elements_with_integrals.h"
#include "modified_volume_constraint_elements_with_integrals.h"
#include "hele_shaw_bubble_problem_with_integrals.h"


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::doc_solution(const std::string& comment)
{
  oomph_info << "Docing step: " << Problem_Parameter::Doc_info.number()
             << std::endl;

  ofstream some_file;
  char filename[100];
  sprintf(filename,
          "%s/soln_%i.dat",
          Problem_Parameter::Doc_info.directory().c_str(),
          Problem_Parameter::Doc_info.number());

  // Number of plot points
  unsigned npts;
  npts = 3;

  // Compute errors and assign to each element for plotting
  double max_err;
  double min_err;
  compute_error_estimate(max_err, min_err);
  //


  some_file.open(filename);
  this->Fluid_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 25, Y = 78, CS=FRAME T = \"Global Step "
            << Problem_Parameter::Doc_info.number() << "  " << comment
            << "\"\n";
  some_file.close();


  // Output boundaries
  sprintf(filename,
          "%s/boundaries_%i.dat",
          Problem_Parameter::Doc_info.directory().c_str(),
          Problem_Parameter::Doc_info.number());
  some_file.open(filename);
  this->Fluid_mesh_pt->output_boundaries(some_file);
  some_file.close();

  // Get max/min area
  double max_area = 0;
  double min_area = 0;
  Fluid_mesh_pt->max_and_min_element_size(max_area, min_area);

  // Write trace file
  Problem_Parameter::Trace_file
    << Problem_Parameter::Doc_info.number() << " " // 1
    << this->time_pt()->time() << " " // 2
    << Fluid_mesh_pt->nelement() << " " // 3
    << max_area << " " // 4
    << min_area << " " // 5
    << get_V() << " " // 6
    << get_P() << " " // 7
    << get_U() << " " // 8
    << get_Q_inv() << " " // 9
    << get_Q() << " " // 10
    << get_Q() * get_U() << " " // 11 Capillary number
    << get_h() << " " // 12
    << get_w() << " " // 13
    << get_alpha() << " " // 14
    << get_asymmetry() << " " // 15
    << get_CoM_Y() << " "; // 16  Bubble centre of mass

  std::cout << "Integral measures" << std::endl;
  for (unsigned i_measure = 0;
       i_measure < Problem_Parameter::n_integral_measures;
       i_measure++)
  {
    std::cout << i_measure << " " << Integral_measures_data_pt->value(i_measure)
              << std::endl;
    Problem_Parameter::Trace_file << Integral_measures_data_pt->value(i_measure)
                                  << " ";
  }
  Problem_Parameter::Trace_file << std::endl;

  // Increment the doc_info number
  Problem_Parameter::Doc_info.number()++;

} // end_of_doc_solution


//==========start_of_main=====================================
/// Driver code for moving bubble problem
//============================================================
int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc, argv);
#endif

  unsigned n_steps = 150;
  double maj_rad;
  double CoM_unsteady;
  double Q_unsteady;
  string output_directory = "";
  bool has_unrecognised_arg = false;

  // REQUIRED
  CommandLineArgs::specify_command_line_flag(
    "-n", &n_steps, "Number of time steps");
  CommandLineArgs::specify_command_line_flag(
    "-r", &maj_rad, "Bubble major radius");
  CommandLineArgs::specify_command_line_flag(
    "-c", &CoM_unsteady, "Bubble centre of mass offset");
  CommandLineArgs::specify_command_line_flag(
    "-q", &Q_unsteady, "Channel flux, Q");

  // OPTIONAL
  CommandLineArgs::specify_command_line_flag(
    "-o",
    &output_directory,
    "Optional: Output directory (e.g. data/bubble_unsteady/ )");

  CommandLineArgs::parse_and_assign(argc, argv, &has_unrecognised_arg);

  if (output_directory != "")
  {
    Problem_Parameter::Doc_info.set_directory(output_directory);
  }
  else
  {
    Problem_Parameter::Doc_info.set_directory("data/bubble_unsteady/");
  }


  // Create generalised Hookean constitutive equations
  Problem_Parameter::Constitutive_law_pt =
    new GeneralisedHookean(&Problem_Parameter::Nu);

  // Open trace file
  Problem_Parameter::Trace_file.open(Problem_Parameter::Doc_info.directory() +
                                     "trace_bubble_test.dat");
  // Increase precision of output
  Problem_Parameter::Trace_file.precision(20);

  // Open norm file
  Problem_Parameter::Norm_file.open(Problem_Parameter::Doc_info.directory() +
                                    "norm.dat");
  Problem_Parameter::OccluHeight_file.open(
    Problem_Parameter::Doc_info.directory() + "Occlusion_Height.dat");
  Problem_Parameter::UpperWall_file.open(
    Problem_Parameter::Doc_info.directory() + "UpperWall_trace.dat");

  Problem_Parameter::Length = 4;
  /// Sets Volume of Initial Condition
  double rad = 0.46;
  double vol = MathematicalConstants::Pi * rad * rad;

  Problem_Parameter::Major_Radius =
    maj_rad; //*0.46//Ellipse1 = 0.5,Ellipse2=0.6,Ellipse3=0.9
  Problem_Parameter::Minor_Radius = vol / (MathematicalConstants::Pi * maj_rad);
  Problem_Parameter::xcenter = 0.0;

  std::cout << " " << std::endl;
  std::cout << "You have chosen ..." << std::endl;
  std::cout << "Maj radius is " << Problem_Parameter::Major_Radius << std::endl;
  std::cout << "Min radius is " << Problem_Parameter::Minor_Radius << std::endl;
  std::cout << "Volume is "
            << MathematicalConstants::Pi * Problem_Parameter::Major_Radius *
                 Problem_Parameter::Minor_Radius
            << std::endl;


  Problem_Parameter::ycenter = CoM_unsteady;
  Problem_Parameter::circpts = 64.0;

  /// The system is currently set up with 12 integral measures (search for
  /// i_measure in hele_shaw_interface_elements_with_integrals.h). But that
  /// number can be changed here
  Problem_Parameter::n_integral_measures = 13;

  BubbleInChannelProblem<MyNewElement> problem;

  maj_rad = problem.get_maj_rad();
  double min_rad = problem.get_min_rad();
  problem.set_w(0.25);
  problem.set_V(MathematicalConstants::Pi * maj_rad * min_rad);
  problem.set_alpha(40);
  problem.set_h(0.024);

  ofstream some_file;
  char filename[100];
  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "initial.dat").c_str());
  some_file.open(filename);
  some_file << Problem_Parameter::Major_Radius << " "
            << Problem_Parameter::Minor_Radius << " "
            << Problem_Parameter::ycenter << " " << problem.get_Q() << "\"\n";
  some_file.close();


  ///   --------------- unsteady solver -----------------
  std::cout << "===========================================" << std::endl;
  std::cout << " INITIAL CONDITION " << std::endl;
  std::cout << "===========================================" << std::endl;

  std::cout << "Bubble pressure : " << problem.get_P() << std::endl;
  std::cout << "Frame speed: " << problem.get_U() << std::endl;
  std::cout << "Bubble volume: " << problem.get_V() << std::endl;
  std::cout << "h is " << problem.get_h() << std::endl;
  std::cout << "Major Radius is " << problem.get_maj_rad() << std::endl;
  std::cout << "Minor Radius is " << problem.get_min_rad() << std::endl;
  std::cout << "X Centre is " << problem.get_xcenter() << std::endl;
  std::cout << "Y Centre is " << problem.get_ycenter() << std::endl;
  std::cout << "No. of points on the circle " << problem.get_circpts()
            << std::endl;
  std::cout << "COM_Y is " << problem.get_CoM_Y() << std::endl;

  problem.doc_solution();

  problem.reset_lagrangian_coordinates();

  double dt = 0.01;
  problem.initialise_dt(dt);
  problem.assign_initial_values_impulsive(dt);

  problem.set_Q(Q_unsteady);

  unsigned max_adapt = 0;
  bool first_step = true;
  for (unsigned i = 0; i < n_steps; i++)
  {
    if (i % 5 == 4)
    {
      max_adapt = 1;
    }
    else
    {
      max_adapt = 0;
    }
    problem.unsteady_newton_solve(dt, max_adapt, first_step);
    std::cout << "===========================================" << std::endl;
    std::cout << " UNSTEADY STEP " << std::endl;
    std::cout << "===========================================" << std::endl;
    first_step = false;
    problem.doc_solution();
    std::cout << "Bubble pressure : " << problem.get_P() << std::endl;
    std::cout << "Frame speed: " << problem.get_U() << std::endl;
    std::cout << "Bubble volume: " << problem.get_V() << std::endl;
    std::cout << "Q is " << problem.get_Q() << std::endl;
    std::cout << "h is " << problem.get_h() << std::endl;
    std::cout << "COM_Y is " << problem.get_CoM_Y() << std::endl;
    problem.reset_lagrangian_coordinates();
  }


#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

  return 0;
} // End of main
