#include <iostream>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"
#include "info_element.h"
#include "my_constraint_elements.h"

#include "inject_air_problem.h"
#include "inject_air_parameters.h"
//#include "custom_hele_shaw_elements_with_integrals.h"

using namespace oomph;
using namespace std;

/// Set parameters
void set_parameters()
{
  /// Set dimensionless parameters
  inject_air::ca_inv = 0.01;
  inject_air::st = 1.0;
  inject_air::alpha = 40;

  // Set bubble target volume
  inject_air::target_bubble_volume_pt = new double;
  *inject_air::target_bubble_volume_pt = -inject_air::initial_volume;

  inject_air::major_radius = 0.05;

  // Create generalised Hookean constitutive equations
  inject_air::nu = 0.3;
  inject_air::constitutive_law_pt =
    new GeneralisedHookean(&inject_air::nu);
}

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

  set_parameters();

  /// Create problem
  RelaxingBubbleProblem<ProjectableTHeleShawPVDElement> problem;

  /// Solve for initial conditions
  problem.solve_for_initial_conditions(doc_info);

  /// Run self tests and Jacobian test
  bool run_self_test = false;
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
  double tF;
  if (CommandLineArgs::command_line_flag_has_been_set(validate_flag_string))
  {
    tF = 5 * dt;
  }
  else
  {
    tF = 1e0;
  }
  problem.iterate_timestepper(dt, tF, doc_info);

  delete inject_air::constitutive_law_pt;
  delete inject_air::target_bubble_volume_pt;

  return 0;
}
