#include <iostream>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

#include "relaxing_bubble_problem.h"
//#include "custom_hele_shaw_elements_with_integrals.h"

using namespace oomph;
using namespace std;

int main(int argc, char* argv[])
{
  cout << "Relaxing bubble demo" << endl;

  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT/");

  relaxing_bubble::total_volume = 1.0;
  relaxing_bubble::target_bubble_volume = MathematicalConstants::Pi *
                                          relaxing_bubble::circular_radius *
                                          relaxing_bubble::circular_radius;
  relaxing_bubble::target_fluid_volume =
    relaxing_bubble::total_volume - relaxing_bubble::target_bubble_volume;

  /// Set parameters
  relaxing_bubble::major_radius = 0.5;
  relaxing_bubble::ca_inv = 10.0;
  relaxing_bubble::st = 1.0;
  relaxing_bubble::alpha = 40.0;
  relaxing_bubble::nu = 0.3;
  // Create generalised Hookean constitutive equations
  relaxing_bubble::constitutive_law_pt =
    new GeneralisedHookean(&relaxing_bubble::nu);

  /// Create problem
  RelaxingBubbleProblem<MyNewElement> problem;

  bool run_self_test = true;
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

  double dt = 5e-3;
  double tF = 2e-1;

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  delete relaxing_bubble::constitutive_law_pt;

  return 0;
}
