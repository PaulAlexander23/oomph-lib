#include <iostream>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"

#include "relaxing_bubble_problem.h"

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
  doc_info.set_directory("RESLT");

  /// Set parameters
  relaxing_bubble::q_inv_pt = new double;
  relaxing_bubble::st_pt = new double;
  relaxing_bubble::alpha_pt = new double;
  *relaxing_bubble::q_inv_pt = 0.01;
  *relaxing_bubble::st_pt = 1;
  *relaxing_bubble::alpha_pt = 0.1;

  /// Create problem
  RelaxingBubbleProblem<ProjectableHeleShawElement<
    PseudoSolidNodeUpdateElement<THeleShawElement<3>, TPVDElement<2, 3>>>>
    problem;

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

  double dt = 0.1;
  double tF = 0.5;

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  return 0;
}
