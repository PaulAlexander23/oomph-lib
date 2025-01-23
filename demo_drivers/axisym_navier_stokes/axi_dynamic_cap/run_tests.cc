#define BOOST_TEST_MODULE Main
#include <boost/test/included/unit_test.hpp>

#include "generic.h"
#include "axisym_navier_stokes.h"
// #include "singular_axisym_dynamic_cap_problem.h"
// #include "parameters.h"
#include "run_tests.h"
#include "full_continuation_problem.h"
#include "my_eigenproblem.h"

using namespace std;
using namespace oomph;

/// ----------------------------------------------------------------------------
/// Matrix

/// Test matrix comparison equal
BOOST_AUTO_TEST_CASE(compare_matrix_same)
{
  const unsigned N = 4;
  DenseDoubleMatrix A(N, N, 0.0);
  DenseDoubleMatrix B(N, N, 0.0);
  A(1, 1) = 1.0;
  B(1, 1) = 1.0;
  B(2, 1) = 1e-10;
  BOOST_TEST(compare_matrices(A, B) == 1);
}

/// Test matrix comparison different
BOOST_AUTO_TEST_CASE(compare_matrix_different)
{
  const unsigned N = 4;
  DenseDoubleMatrix A(N, N, 0.0);
  DenseDoubleMatrix B(N, N, 0.0);
  A(2, 1) = 1.0;
  B(2, 1) = 1e-10;
  BOOST_TEST(compare_matrices(A, B) == 0);
}

/// ----------------------------------------------------------------------------
/// Problem

/// Test temp
BOOST_AUTO_TEST_CASE(temp)
{
  Temp temp;
}

/// Test MyEigenproblem
BOOST_AUTO_TEST_CASE(my_eigenproblem_creation)
{
  MyEigenproblem problem;
}

/// Test SingularAxisymDynamicCapProblem
BOOST_AUTO_TEST_CASE(acute_axisym_problem)
{
  Params parameters;
  parameters.contact_angle = 60.0 * MathematicalConstants::Pi / 180.0;
  AXISYM_PROBLEM problem(&parameters);
  problem.setup();
  //BOOST_TEST(problem.debug_jacobian());
  problem.newton_solve();
  // If the problem solve completes, then it is a success.
}

/// Test SingularAxisymDynamicCapProblem
BOOST_AUTO_TEST_CASE(singular_axisym_dynamic_cap_problem_creation)
{
  Params parameters;
  AXISYM_PROBLEM problem(&parameters);
  problem.setup();
  //BOOST_TEST(problem.debug_jacobian());
  problem.newton_solve();
  // If the problem solve completes, then it is a success.
}

/// ----------------------------------------------------------------------------
/// Continuation

/// Test adaptive step size
BOOST_AUTO_TEST_CASE(height_continuation_problem_adapt_size)
{
  HeightContinuationProblem problem;
  // This is the default value
  const unsigned desired_newton_iterations_ds = 5;
  problem.set_nnewton_iter_taken(1);
  double actual = problem.adapt_step_size(1.0);
  BOOST_TEST(actual == 1.5);

  problem.set_nnewton_iter_taken(desired_newton_iterations_ds + 1);
  actual = problem.adapt_step_size(1.0);
  BOOST_TEST(actual == 2.0 / 3.0);

  problem.set_nnewton_iter_taken(desired_newton_iterations_ds);
  actual = problem.adapt_step_size(1.0);
  BOOST_TEST(actual == 1.0);
}

/// Test the height doc solution
BOOST_AUTO_TEST_CASE(height_continuation_problem_doc_solution)
{
  HeightContinuationProblem problem;
  problem.doc_solution();
}

/// Test straight-forward continuation
BOOST_AUTO_TEST_CASE(continuation)
{
  Params parameters;
  SingularAxisymDynamicCapProblem<
    SingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>
    problem(&parameters);
  problem.setup();

  double ds = 0.1;
  ds = problem.continuation_step_solve(
    parameters.reynolds_inverse_froude_number_pt, ds);

  BOOST_TEST(ds > 0);
}

/// Full continuation problem creation
// BOOST_AUTO_TEST_CASE(full_continuation_problem_reynolds_inverse_froude_number)
//{
//  Params parameters;
//  FullContinuationProblem<BASE_ELEMENT, TIMESTEPPER> problem(parameters);
//  problem.setup();
//  problem.set_continuation_parameter(parameters.reynolds_inverse_froude_number);
//  BOOST_TEST(problem.debug_jacobian());
//  // problem.newton_solve();
//  //  If the problem solve completes, then it is a success.
//}
