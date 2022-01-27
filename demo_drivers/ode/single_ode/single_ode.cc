#include <iostream>

#include "generic.h"

using namespace std;
using namespace oomph;

#include "ode_problem_inbuilt.h"

int main(int argc, char* argv[])
{
  CommandLineArgs::setup(argc, argv);

  string output_directory = "RESLT";
  CommandLineArgs::specify_command_line_flag(
    "-o",
    &output_directory,
    "Optional: Output directory without the trailing slash (e.g. RESLT )");

  bool has_unrecognised_arg = false;
  CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);

  DocInfo doc_info(output_directory);

  ODEProblem my_ode_problem;

  my_ode_problem.set_initial_condition();
  my_ode_problem.doc_solution(doc_info);

  double t_step = 0.1;
  double t_final = 2.0;
  my_ode_problem.iterate_timestepper(t_step, t_final, doc_info);

  return 0;
}
