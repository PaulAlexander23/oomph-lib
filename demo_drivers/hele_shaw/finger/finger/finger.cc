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
  const double Q = 0.05;
  const double circular_radius = 0.0;
  const double V =
    MathematicalConstants::Pi * circular_radius * circular_radius;
  const double h = 0.024;

  double Ca;
  double new_V;
  double new_h;
  convert_parameters_from_q_nd_to_capillary_nd(Q, V, h, Ca, new_V, new_h);

  finger::alpha = 40.0;
  finger::ca_inv = 1.0 / Ca;
  finger::st = 1.0;
  finger::nu = 0.3;

  finger::major_radius = 0.25;
  finger::bubble_initial_centre_y = 0.5 + 0.005;
  // finger::target_bubble_volume = 0.01; // new_volume;

  finger::perturbation_amplitude = new_h;
  finger::perturbation_rms_width = 0.25;

  // Create generalised Hookean constitutive equations
  finger::constitutive_law_pt = new GeneralisedHookean(&finger::nu);

  // IntegralProblem<QIntegralElement<3>>
  // integral_problem(finger::channel_depth); integral_problem.newton_solve();
  // integral_problem.doc_solution(doc_info);
  // finger::total_volume = integral_problem.result();

  finger::target_bubble_volume =
    MathematicalConstants::Pi * pow(finger::finger_width, 2.0) / 2;
  finger::total_volume = 1.0;

  finger::target_fluid_volume =
    (finger::total_volume) - (finger::target_bubble_volume);

  finger::total_flux = 1.0;

  /// Print the parameters
  finger::print_parameters();

  /// Create problem
  FingerProblem<MyNewElementWithIntegral> problem(doc_info);

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

  double dt = 1e-4;
  double tF = 1e-2;

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  delete finger::constitutive_law_pt;

  return 0;
}
