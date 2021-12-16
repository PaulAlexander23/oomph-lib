#include <iostream>
#include <chrono>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

#include "moving_bubble_problem.h"
#include "integral_problem.h"

using namespace oomph;
using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])
{
  auto start_time = high_resolution_clock::now();
  cout << "Moving bubble demo with moving frame" << endl;

  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT/");

  /// Set parameters
  const double major_radius = 0.5;
  const double width = 0.25; // 2;
  const double circular_radius = 0.46;
  const double volume =
    MathematicalConstants::Pi * circular_radius * circular_radius;
  double Q = 0.05;

  CommandLineArgs::specify_command_line_flag("-q", &Q, "Channel flux, Q");

  double length_ratio;
  double pressure_ratio;
  double velocity_ratio;
  double time_ratio;
  double new_r;
  double new_width;
  double new_volume;
  double new_Ca;
  moving_bubble::alpha = 40.0;

  convert_to_capillary_nondimensionalisation(major_radius,
                                             width,
                                             volume,
                                             Q,
                                             moving_bubble::alpha,
                                             length_ratio,
                                             pressure_ratio,
                                             velocity_ratio,
                                             time_ratio,
                                             new_r,
                                             new_width,
                                             new_volume,
                                             new_Ca);
  moving_bubble::major_radius = new_r;
  moving_bubble::ca_inv = 0.5;//1.0 / new_Ca;
  moving_bubble::st = 1.0;
  moving_bubble::nu = 0.3;

  moving_bubble::bubble_initial_centre_y = 0.5 * (1 + 0.01);

  moving_bubble::perturbation_amplitude = 0.024;
  moving_bubble::perturbation_rms_width = new_width;

  moving_bubble::target_bubble_volume = new_volume;
  // Create generalised Hookean constitutive equations
  moving_bubble::constitutive_law_pt =
    new GeneralisedHookean(&moving_bubble::nu);

  moving_bubble::global_frame_travel_pt = new double;
  *moving_bubble::global_frame_travel_pt = 0.0;

  IntegralProblem<QIntegralElement<3>> integral_problem(
    moving_bubble::channel_depth);
  integral_problem.newton_solve();
  moving_bubble::total_volume = integral_problem.result();
  integral_problem.doc_solution(doc_info);
  cout << "Total volume: " << moving_bubble::total_volume << endl;


  moving_bubble::target_fluid_volume =
    (moving_bubble::total_volume) - (moving_bubble::target_bubble_volume);

  /// Print the parameters
  moving_bubble::print_parameters();

  /// Create problem
  // RelaxingBubbleProblem<ProjectableHeleShawElementWithSolidFaces<3>> problem;
  // RelaxingBubbleProblem<HeleShawWithErrorElement> problem;
  RelaxingBubbleProblem<MyNewElementWithIntegral> problem(doc_info);

  bool run_self_test = false;
  if (run_self_test)
  {
    if (problem.self_test())
    {
      throw OomphLibError(
        "Self test failed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  /// Solve for initial conditions
  problem.solve_for_initial_conditions(doc_info);

  double dt = 3e-2;
  double tF = 3e0;

  // Problem* problem_pt = new RelaxingBubbleProblem<MyNewElementWithIntegral>;
  // DoubleVector result;
  // FD_LU fd_lu;
  // fd_lu.solve(problem_pt, result);
  // cout << "solved" << endl;
  // result.output(cout);

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  delete moving_bubble::constitutive_law_pt;

  auto end_time = high_resolution_clock::now();
  duration<double> diff = duration_cast<minutes>(end_time - start_time);
  cout << "Total real-time duration: " << diff.count() << endl;

  return 0;
}
