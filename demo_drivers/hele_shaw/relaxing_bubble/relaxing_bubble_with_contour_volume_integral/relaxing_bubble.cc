#include <iostream>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"
#include "info_element.h"
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
  relaxing_bubble::ca_inv = 100;
  relaxing_bubble::st = 1.0;
  relaxing_bubble::alpha = 40;

  // Set bubble target volume
  relaxing_bubble::target_bubble_volume =
    -MathematicalConstants::Pi * pow(relaxing_bubble::circular_radius, 2.0);

  relaxing_bubble::major_radius = 0.25;

  // Create generalised Hookean constitutive equations
  relaxing_bubble::nu = 0.3;
  relaxing_bubble::constitutive_law_pt =
    new GeneralisedHookean(&relaxing_bubble::nu);
}

int main(int argc, char* argv[])
{
  cout << "Relaxing bubble demo" << endl;

  /// Setup and store command line arguments
  string validate_flag_string = "--validate";
  string self_test_flag_string = "--self_test";
  bool has_unrecognised_arg = false;
  CommandLineArgs::setup(argc, argv);
  CommandLineArgs::specify_command_line_flag(validate_flag_string,
                                             "Optional: Run with self tests.");
  CommandLineArgs::specify_command_line_flag(
    self_test_flag_string, "Optional: Run with validation parameters.");
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
  if (CommandLineArgs::command_line_flag_has_been_set(self_test_flag_string))
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
  double dt = 4e-2;
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

  delete relaxing_bubble::constitutive_law_pt;

  return 0;
}
