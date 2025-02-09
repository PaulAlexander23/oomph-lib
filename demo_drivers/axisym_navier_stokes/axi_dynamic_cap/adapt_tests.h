#include <boost/test/included/unit_test.hpp>

#include "generic.h"
#include "axisym_navier_stokes.h"
#include "height_continuation_problem.h"
// #include "singular_axisym_dynamic_cap_problem.h"
// #include "parameters.h"
#include "run_tests.h"
#include "full_continuation_problem.h"
#include "my_eigenproblem.h"

BOOST_AUTO_TEST_CASE(initial_mesh_number_of_points)
{
  Params parameters;
  parameters.polyline_refinement_tolerence = 4e-3;
  SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> problem(
    &parameters);
  problem.newton_solve();
  double max_free_surface_error = problem.max_free_surface_error();
  BOOST_TEST(max_free_surface_error == 0.0);
}
