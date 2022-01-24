#include <iostream>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"


#include "inject_air_problem.h"
#include "integral_problem.h"
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

  /// Set parameters
  parameters::ca_inv_pt = new double;
  parameters::st_pt = new double;
  parameters::alpha_pt = new double;
  parameters::nu_pt = new double;
  parameters::target_bubble_volume_pt = new double;
  parameters::target_fluid_volume_pt = new double;
  parameters::total_volume_pt = new double;
  *parameters::ca_inv_pt = 0.141376;
  *parameters::st_pt = 1.0;
  *parameters::alpha_pt = 40.0;
  *parameters::nu_pt = 0.3;
  *parameters::target_bubble_volume_pt = parameters::initial_volume;
  // Create generalised Hookean constitutive equations
  parameters::constitutive_law_pt = new GeneralisedHookean(parameters::nu_pt);

  IntegralProblem<QIntegralElement<3>> integral_problem(
    parameters::channel_depth);
  integral_problem.newton_solve();
  *parameters::total_volume_pt = integral_problem.result();
  integral_problem.doc_solution(doc_info);

  *parameters::total_volume_pt = MathematicalConstants::Pi * pow(1.0, 2.0);

  *parameters::target_fluid_volume_pt =
    (*parameters::total_volume_pt) - (*parameters::target_bubble_volume_pt);

  parameters::injection_rate = 1.0;

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

  delete parameters::ca_inv_pt;
  delete parameters::st_pt;
  delete parameters::alpha_pt;
  delete parameters::nu_pt;
  delete parameters::target_bubble_volume_pt;
  delete parameters::target_fluid_volume_pt;
  delete parameters::total_volume_pt;
  delete parameters::constitutive_law_pt;

  return 0;
}
