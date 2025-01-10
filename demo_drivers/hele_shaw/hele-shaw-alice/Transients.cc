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


// The equations
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;

#include "hele_shaw_interface_elements_2.h"
#include "Thele_shaw_elements.h"
#include "hele_shaw_flux_elements.h"
#include "custom_hele_shaw_elements_2.h"
#include "Polygon_manipulation_3.h"
#include "modified_volume_constraint_elements_2.h"
#include "hele_shaw_multi_drops_problem_3.h"

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::write_output_file()
{
  char filename[100];
  sprintf(filename,
          "%s/ordered_%s_%i.dat",
          Problem_Parameter::Doc_info.directory().c_str(),
          Problem_Parameter::Output_label.c_str(),
          Problem_Parameter::Doc_info.number());


  ofstream some_file;
  some_file.open(filename);

  some_file << Problem_Parameter::N_Bubble << std::endl;
  some_file << Problem_Parameter::Doc_info.number() << " "
            << this->time_pt()->time() << std::endl;
  some_file << get_Q() << " " << get_U() << " " << get_h() << " " << get_w()
            << " " << get_alpha() << " " << std::endl;
  some_file << get_drop_viscosity() << std::endl;
  some_file << Problem_Parameter::branch_no << std::endl;


  output_ordered_polygon(some_file);
  some_file.close();
  return;


  unsigned N_0 = Problem_Parameter::Ordered_bubbles.size();
  check_topology();
  std::cout << "Ordered bubbles " << N_0 << std::endl;

  for (unsigned i_bubble = 0; i_bubble < N_0; i_bubble++)
  {
    unsigned i_vertex = Problem_Parameter::Ordered_bubbles[i_bubble].size();

    for (unsigned i = 0; i < i_vertex; i++)
    {
      double x_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][0];
      double y_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][1];
      some_file << x_i << " " << y_i << " " << i_bubble << " " << i_vertex
                << std::endl;
    }
  }


  some_file.close();
}


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::doc_solution(const std::string& comment)
{
  oomph_info << "Docing step: " << Problem_Parameter::Doc_info.number()
             << std::endl;

  char output_filename[100];
  sprintf(output_filename,
          "%s/output_%s_%i.dat",
          Problem_Parameter::Doc_info.directory().c_str(),
          Problem_Parameter::Output_label.c_str(),
          Problem_Parameter::Doc_info.number());


  write_output_file();


  unsigned npts;
  npts = 3;

  double max_err;
  double min_err;
  compute_error_estimate(max_err, min_err);

  ofstream some_file;
  char filename[100];
  sprintf(filename,
          "%s/soln_%s_%i.dat",
          Problem_Parameter::Doc_info.directory().c_str(),
          Problem_Parameter::Output_label.c_str(),
          Problem_Parameter::Doc_info.number());


  some_file.open(filename);
  this->Fluid_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 25, Y = 78, CS=FRAME T = \"Global Step "
            << Problem_Parameter::Doc_info.number() << "  " << comment
            << "\"\n";
  some_file.close();

  sprintf(filename,
          "%s/boundaries_%s_%i.dat",
          Problem_Parameter::Doc_info.directory().c_str(),
          Problem_Parameter::Output_label.c_str(),
          Problem_Parameter::Doc_info.number());
  some_file.open(filename);
  this->Fluid_mesh_pt->output_boundaries(some_file);
  some_file.close();

  double max_area = 0;
  double min_area = 0;
  Fluid_mesh_pt->max_and_min_element_size(max_area, min_area);

  Problem_Parameter::Trace_file
    << Problem_Parameter::Doc_info.number() << " " // 1
    << this->time_pt()->time() + Problem_Parameter::Time_shift << " " // 2
    << this->time_pt()->time() << " " // 3
    << Fluid_mesh_pt->nelement() << " " // 4
    << max_area << " " // 5
    << min_area << " " // 6
    << get_U() << " " // 7
    << get_Q_inv() << " " // 8
    << get_Q() << " " // 9
    << get_h() << " " // 10
    << get_w() << " " // 11
    << get_alpha() << " " // 12
    << get_asymmetry() << " " // 13
    << Problem_Parameter::Min_distance << " " // 14
    << Problem_Parameter::branch_no << " " // 15
    << get_drop_viscosity() << " " // 16
    << Problem_Parameter::Valid_output_point << " "; // 17

  Problem_Parameter::Trace_file << std::endl;

  // Increment the doc_info number
  Problem_Parameter::Doc_info.number()++;
  std::cout << "Finished docing " << std::endl;
} // end_of_doc_solution


//==========start_of_main=====================================
/// Driver code for moving bubble problem
//============================================================
int main(int argc, char** argv)
{
#ifdef OOMPH_HAS_MPI
  //   MPI_Helpers::init(argc,argv);
#endif


  // Store command line arguments
  CommandLineArgs::setup(argc, argv);
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  // Create generalised Hookean constitutive equations
  Problem_Parameter::Constitutive_law_pt =
    new GeneralisedHookean(&Problem_Parameter::Nu);

  Problem_Parameter::Output_label = "Obstacle_transient_visc_0p5_h_0p05";
  Problem_Parameter::Trace_file.open("RESLT_map/trace_" +
                                     Problem_Parameter::Output_label + ".dat");
  Problem_Parameter::Trace_file.precision(20);
  Problem_Parameter::Norm_file.open("RESLT_map/norm.dat");
  Problem_Parameter::Doc_info.set_directory("RESLT_map");

  Problem_Parameter::Channel_start = -3;
  Problem_Parameter::Channel_end = 3;

  double R1 = 0, R2 = 0;

  double X1 = 0, X2 = 0;
  double Y1 = 0, Y2 = 0;

  /// These are flags for the constructor.
  Problem_Parameter::Read_from_file = false;
  Problem_Parameter::Reload_from_vector = true;

  Problem_Parameter::N_Bubble = 1;

  X1 = 0;
  R1 = 0.7;
  Y1 = 0.5;

  double Q = 0.02;
  double e_vertical = 0.5;
  double drop_viscosity = 0.5;
  double obstacle_height = 0.05;


  /// Begin setup for constructor.
  Problem_Parameter::All_the_bubbles.resize(0);
  Problem_Parameter::Bubble_volume_vector.resize(0);
  Problem_Parameter::All_the_bubbles.resize(Problem_Parameter::N_Bubble);
  Problem_Parameter::Bubble_volume_vector.resize(Problem_Parameter::N_Bubble);

  for (unsigned i_bubble = 0; i_bubble < Problem_Parameter::N_Bubble;
       i_bubble++)
  {
    Vector<double> new_point(2, 0.0);
    unsigned n_points = 50;

    double x_center = 0;
    double y_center = 0;
    double radius = 0.0; /// This default value would lead to an error.
    Problem_Parameter::e_vertical = e_vertical;
    if (i_bubble == 0)
    {
      x_center = X1;
      radius = R1;
      y_center = Y1;
    }
    if (i_bubble == 1)
    {
      x_center = X2;
      radius = R2;
      y_center = Y2;
    }

    for (unsigned i_point = 0; i_point < n_points; i_point++)
    {
      double theta = double(i_point) / double(n_points - 1) * 2.0 * M_PI + M_PI;
      double x = x_center + std::cos(theta) * radius / e_vertical;
      double y = y_center + std::sin(theta) * radius * e_vertical;
      new_point[0] = x;
      new_point[1] = y;
      Problem_Parameter::All_the_bubbles[i_bubble].push_back(new_point);
    }
    Problem_Parameter::Bubble_volume_vector[i_bubble] = -M_PI * radius * radius;
  }


  /// Each change in topology requires the construction of a new problem.
  /// So in a loop, we construct a problem and timestep until if/when there is a
  /// change in topology.
  for (unsigned topology_trial = 0; topology_trial < 1; topology_trial++)
  {
    std::cout << "Constructing a new problem" << std::endl;

    Problem_Parameter::Topology_change_needed = false;

    std::cout << "Should have " << Problem_Parameter::N_Bubble << " bubbles "
              << std::endl;

    BubbleInChannelProblem<MyNewElement> problem;

    problem.set_w(0.25);
    problem.set_alpha(40);
    problem.set_h(obstacle_height);
    problem.set_U(0);
    problem.set_drop_viscosity(drop_viscosity);
    problem.set_Q(Q);
    problem.set_G(1);

    problem.unpin_U();
    problem.redo_equation_numbering();
    problem.snap_bubble_to_ellipse(0,
                                   0,
                                   Y1,
                                   R1 / Problem_Parameter::e_vertical,
                                   R1 * Problem_Parameter::e_vertical);
    problem.reset_lagrangian_coordinates();

    double dt = 0.05 * (1 + drop_viscosity);

    problem.initialise_dt(dt);
    problem.assign_initial_values_impulsive(dt);
    problem.doc_solution();
    for (unsigned i = 0; i < 500; i++)
    {
      unsigned max_adapt = 0;
      if (i % 3 == 2)
      {
        max_adapt = 1;
      }
      problem.unsteady_newton_solve(dt, max_adapt, false);
      problem.doc_solution();
      std::cout << "U: " << problem.get_U() << std::endl;
      std::cout << "Time: " << problem.get_time() << std::endl;
      problem.reset_lagrangian_coordinates();
    }
  }


  // Kill const eqn
  delete Problem_Parameter::Constitutive_law_pt;


#ifdef OOMPH_HAS_MPI
  //   MPI_Helpers::finalize();
#endif

  return 0;
} // End of main
