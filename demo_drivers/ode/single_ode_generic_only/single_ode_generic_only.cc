#include <iostream>

#include "generic.h"

#include "ode_elements.h"

using namespace std;
using namespace oomph;

#include "ode_problem_my_element.h"

int main(int argc, char *argv[]) {
  cout << "driver_ode_my_element.cc" << endl;
  CommandLineArgs::setup(argc, argv);

  DocInfo doc_info("data");
  doc_info.enable_error_if_directory_does_not_exist();

  ODEProblem my_ode_problem;

  my_ode_problem.set_initial_condition();
  my_ode_problem.doc_solution(doc_info);

  double t_step = 0.5;
  double t_final = 2.0;
  my_ode_problem.iterate_timestepper(t_step, t_final, doc_info);

  return 0;
}
