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

  /// Setup and store command line arguments
  string validate_flag_string = "--validate";
  bool has_unrecognised_arg = false;
  CommandLineArgs::setup(argc, argv);
  CommandLineArgs::specify_command_line_flag(validate_flag_string,
                                             "Optional: Run with self tests.");
  CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT/");

  /// Set parameters
  parameters::target_bubble_volume_pt = new double;
  parameters::target_fluid_volume_pt = new double;
  parameters::total_volume_pt = new double;
  parameters::ca_inv = 0.01;
  parameters::st = 1.0;
  parameters::alpha = 40.0;
  parameters::nu = 0.3;
  *parameters::target_bubble_volume_pt = parameters::initial_volume;
  // Create generalised Hookean constitutive equations
  parameters::constitutive_law_pt = new GeneralisedHookean(&parameters::nu);

  IntegralProblem<QIntegralElement<3>> integral_problem(
    parameters::channel_depth);
  integral_problem.newton_solve();
  *parameters::total_volume_pt = integral_problem.result();
  integral_problem.doc_solution(doc_info);

  *parameters::total_volume_pt = MathematicalConstants::Pi * pow(1.0, 2.0);

  *parameters::target_fluid_volume_pt =
    (*parameters::total_volume_pt) - (*parameters::target_bubble_volume_pt);

  parameters::injection_rate = 0.1;

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

  double dt = 5e-3;
  double tF;
  if (CommandLineArgs::command_line_flag_has_been_set(validate_flag_string))
  {
    tF = 5 * dt;
  }
  else
  {
    tF = 1e0;
  }

  // Problem* problem_pt = new RelaxingBubbleProblem<MyNewElementWithIntegral>;
  // DoubleVector result;
  // FD_LU fd_lu;
  // fd_lu.solve(problem_pt, result);
  // cout << "solved" << endl;
  // result.output(cout);

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  delete parameters::target_bubble_volume_pt;
  delete parameters::target_fluid_volume_pt;
  delete parameters::total_volume_pt;
  delete parameters::constitutive_law_pt;

  return 0;
}
