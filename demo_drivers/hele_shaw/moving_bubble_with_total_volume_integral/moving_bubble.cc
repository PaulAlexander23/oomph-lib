#include <iostream>

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

int main(int argc, char* argv[])
{
  cout << "Moving bubble demo" << endl;

  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT/");

  /// Set parameters
  const double major_radius = 0.5;
  const double width = 2;
  const double circular_radius = 0.46;
  const double volume =
    MathematicalConstants::Pi * circular_radius * circular_radius;
  const double Q = 0.05;

  double length_ratio;
  double pressure_ratio;
  double velocity_ratio;
  double time_ratio;
  double new_r;
  double new_width;
  double new_volume;
  double new_Ca;

  convert_to_capillary_nondimensionalisation(major_radius,
                                             width,
                                             volume,
                                             Q,
                                             length_ratio,
                                             pressure_ratio,
                                             velocity_ratio,
                                             time_ratio,
                                             new_r,
                                             new_width,
                                             new_volume,
                                             new_Ca);
  moving_bubble::major_radius = new_r;
  moving_bubble::ca_inv = 1.0 / new_Ca;
  moving_bubble::ca_inv = 20.0;
  moving_bubble::st = 1.0;
  moving_bubble::alpha = 40.0;
  moving_bubble::nu = 0.3;

  moving_bubble::bubble_initial_centre_y = 0.5 + 0.005;

  moving_bubble::perturbation_amplitude = 0.1;//0.024;
  moving_bubble::perturbation_rms_width = new_width;

  moving_bubble::target_bubble_volume = new_volume;
  // Create generalised Hookean constitutive equations
  moving_bubble::constitutive_law_pt =
    new GeneralisedHookean(&moving_bubble::nu);

  /// Test with no flux
  total_flux = 0.0;

  IntegralProblem<QIntegralElement<3>> integral_problem(
    moving_bubble::channel_depth);
  integral_problem.newton_solve();
  moving_bubble::total_volume = integral_problem.result();
  integral_problem.doc_solution(doc_info);


  moving_bubble::target_fluid_volume =
    (moving_bubble::total_volume) - (moving_bubble::target_bubble_volume);

  /// Print the parameters
  moving_bubble::print_parameters();

  /// Create problem
  // RelaxingBubbleProblem<ProjectableHeleShawElementWithSolidFaces<3>> problem;
  // RelaxingBubbleProblem<HeleShawWithErrorElement> problem;
  RelaxingBubbleProblem<MyNewElementWithIntegral> problem;

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

  double dt = 1e-2;
  double tF = 1e0;

  // Problem* problem_pt = new RelaxingBubbleProblem<MyNewElementWithIntegral>;
  // DoubleVector result;
  // FD_LU fd_lu;
  // fd_lu.solve(problem_pt, result);
  // cout << "solved" << endl;
  // result.output(cout);

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  delete moving_bubble::constitutive_law_pt;

  return 0;
}
