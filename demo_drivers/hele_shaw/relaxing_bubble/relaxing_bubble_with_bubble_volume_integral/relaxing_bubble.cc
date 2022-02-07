#include <iostream>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"
#include "info_element.h"
#include "projectable_integral_element.h"
#include "Tintegral_elements.h"
#include "integral_elements.h"
#include "my_constraint_elements.h"

#include "relaxing_bubble_problem.h"
#include "relaxing_bubble_parameters.h"
//#include "custom_hele_shaw_elements_with_integrals.h"

using namespace oomph;
using namespace std;

/// Set parameters
void set_parameters()
{
  /// Set dimensionless parameters
  relaxing_bubble::ca_inv = 0.01;
  relaxing_bubble::st = 1.0;
  relaxing_bubble::alpha = 40;

  // Set bubble target volume
  relaxing_bubble::target_bubble_volume_pt =
    new double(4.0 * atan(1.0) * pow(0.3, 2.0));

  // Create generalised Hookean constitutive equations
  relaxing_bubble::nu = 0.3;
  relaxing_bubble::constitutive_law_pt =
    new GeneralisedHookean(&relaxing_bubble::nu);
}

int main(int argc, char* argv[])
{
  cout << "Relaxing bubble demo" << endl;

  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT/");

  set_parameters();

  /// Create problem
  RelaxingBubbleProblem<
    MyNewElement,
    ProjectableIntegralElement<
      PseudoSolidNodeUpdateElement<TIntegralElement<3>, TPVDElement<2, 3>>>>
    problem;

  /// Solve for initial conditions
  problem.solve_for_initial_conditions(doc_info);

  /// Run self tests and Jacobian test
  bool run_self_test = true;
  if (run_self_test)
  {
    bool self_test_failed = problem.self_test();
    if (self_test_failed)
    {
      throw OomphLibError(
        "Self test failed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    problem.debug_jacobian();
  }

  /// Iterate the timestepper using the fixed time step until the final time
  double dt = 5e-3;
  double tF = 2e-1;
  problem.iterate_timestepper(dt, tF, doc_info);

  delete relaxing_bubble::constitutive_law_pt;

  return 0;
}
