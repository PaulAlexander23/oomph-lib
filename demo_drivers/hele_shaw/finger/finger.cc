#include <iostream>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

#include "finger_problem.h"
#include "integral_problem.h"

using namespace oomph;
using namespace std;

int main(int argc, char* argv[])
{
  cout << "Finger demo" << endl;

  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT/");

  /// Set parameters
  const double major_radius = 0.5;
  const double width = 0.25;
  const double circular_radius = 0.0;
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
  finger::major_radius = new_r;
  finger::ca_inv = 1.0 / new_Ca;
  finger::st = 1.0;
  finger::alpha = 40.0;
  finger::nu = 0.3;

  finger::bubble_initial_centre_y = 0.5 + 0.005;

  finger::perturbation_amplitude = 0.0;
  finger::perturbation_rms_width = new_width;

  finger::target_bubble_volume = new_volume;
  // Create generalised Hookean constitutive equations
  finger::constitutive_law_pt = new GeneralisedHookean(&finger::nu);

  //IntegralProblem<QIntegralElement<3>> integral_problem(finger::channel_depth);
  //integral_problem.newton_solve();
  //integral_problem.doc_solution(doc_info);
  //finger::total_volume = integral_problem.result();

  finger::total_volume = 1.0;

  finger::target_fluid_volume =
    (finger::total_volume) - (finger::target_bubble_volume);

  finger::total_flux = 0.0;

  /// Print the parameters
  finger::print_parameters();

  /// Create problem
  FingerProblem<MyNewElementWithIntegral> problem;

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

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  delete finger::constitutive_law_pt;

  return 0;
}
