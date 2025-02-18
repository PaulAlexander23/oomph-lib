#define BOOST_TEST_MODULE axisym_cap_test_module
#include <boost/test/included/unit_test.hpp>

#include "run_tests.h"

using namespace std;
using namespace oomph;


// A test to show the artefact
BOOST_AUTO_TEST_CASE(show_artefact_mode_0)
{
  Params parameters;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  base_problem.steady_newton_solve();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();
  perturbed_problem.doc_solution();

  // Steady newton solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  if (parameters.azimuthal_mode_number != 0)
  {
    perturbed_problem.make_unsteady();
  }
  perturbed_problem.pin_horizontal_mesh_deformation();


  // Eigensolve
  Vector<std::complex<double>> eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(eigenvalue[0].real() - (-1.5677282757727009)) < 1e-6);
}

BOOST_AUTO_TEST_CASE(show_artefact_mode_1)
{
  Params parameters;
  parameters.azimuthal_mode_number = 1;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  base_problem.steady_newton_solve();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();
  perturbed_problem.doc_solution();

  // Steady newton solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  perturbed_problem.make_unsteady();
  perturbed_problem.pin_horizontal_mesh_deformation();


  // Eigensolve
  Vector<std::complex<double>> eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(eigenvalue[0].real() - (-0.59043964949838124)) < 1e-6);
}

BOOST_AUTO_TEST_CASE(show_artefact_mode_2)
{
  Params parameters;
  parameters.azimuthal_mode_number = 2;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  base_problem.steady_newton_solve();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();
  perturbed_problem.doc_solution();

  // Steady newton solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  perturbed_problem.make_unsteady();
  perturbed_problem.pin_horizontal_mesh_deformation();


  // Eigensolve
  Vector<std::complex<double>> eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(eigenvalue[0].real() - (-0.59043964949838124)) < 1e-6);
}


// Mode 0
// No displacement
// BOOST_AUTO_TEST_CASE(mode_0_no_displacement)
// {
//   // Create the parameters
//   Params parameters;
//   parameters.azimuthal_mode_number = 0;
//   parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;
//   //*parameters.wall_velocity_pt = 0.01;
// 
//   // Create the base problem
//   BASE_PROBLEM base_problem(&parameters);
//   base_problem.steady_newton_solve();
//   base_problem.reset_lagrange();
//   base_problem.assign_initial_values_impulsive();
// 
//   // Create the linear problem
//   PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
//                                       base_problem.free_surface_mesh_pt(),
//                                       base_problem.slip_surface_mesh_pt(),
//                                       &parameters);
// 
//   perturbed_problem.pin_horizontal_mesh_deformation();
//   perturbed_problem.pin_vertical_mesh_deformation();
//   perturbed_problem.set_constant_lagrange_free_surface_boundary_condition(0.0,
//                                                                           0.0);
//   perturbed_problem.assign_initial_values_impulsive();
//   perturbed_problem.doc_solution();
// 
//   // Steady newton solve
//   perturbed_problem.steady_newton_solve();
//   perturbed_problem.doc_solution();
// 
//   perturbed_problem.pin_horizontal_mesh_deformation();
// 
//   // Eigensolve
//   Vector<std::complex<double>> eigenvalue =
//     perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
//   // BOOST_TEST(abs(eigenvalue[0].real() - (8.046912588166538e-05)) < 1e-6);
// }

// Augmented linear problem
BOOST_AUTO_TEST_CASE(augmented_linear_problem)
{
  // Create the parameters
  Params parameters;
  parameters.azimuthal_mode_number = 1;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;
  //*parameters.wall_velocity_pt = 0.01;

  // Create the base problem
  typedef SolidSingularAxisymNavierStokesElement<
    ProjectableAxisymmetricTTaylorHoodPVDElement>
    BASE_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> base_problem(
    &parameters);
  base_problem.steady_newton_solve();
  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  parameters.azimuthal_mode_number = 0;
  unsigned local_doc_number = 0;

  // Create the linear problem
  typedef SingularOverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  SingularPerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                             PERTURBED_ELEMENT,
                                             TIMESTEPPER>
    perturbed_problem(base_problem.bulk_mesh_pt(),
                      base_problem.free_surface_mesh_pt(),
                      base_problem.slip_surface_mesh_pt(),
                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.doc_solution();

  // Steady newton solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // perturbed_problem.perturb_vertical_velocity();
  // perturbed_problem.assign_initial_values_impulsive();
  // perturbed_problem.doc_solution();

  //// Time step (docs within the timestepper)
  // perturbed_problem.timestep(0.01, 0.01);

  // perturbed_problem.make_unsteady();
  perturbed_problem.pin_horizontal_mesh_deformation();


  // Eigensolve
  Vector<std::complex<double>> eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(eigenvalue[0].real() - (-1.5677282757561151)) < 1e-6);

  local_doc_number = perturbed_problem.doc_info().number();

  // debug_jacobian(&perturbed_problem);

  parameters.azimuthal_mode_number = 1;

  // Create the linear problem
  SingularPerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                             PERTURBED_ELEMENT,
                                             TIMESTEPPER>
    perturbed_problem0(base_problem.bulk_mesh_pt(),
                       base_problem.free_surface_mesh_pt(),
                       base_problem.slip_surface_mesh_pt(),
                       &parameters);
  perturbed_problem0.assign_initial_values_impulsive();
  perturbed_problem0.doc_info().number() = local_doc_number;
  perturbed_problem0.steady_newton_solve();
  perturbed_problem0.make_unsteady();
  perturbed_problem0.pin_horizontal_mesh_deformation();
  eigenvalue =
    perturbed_problem0.solve_and_document_n_most_unstable_eigensolutions(1);
  BOOST_TEST(abs(eigenvalue[0].real() - (-0.59043964949883476)) < 1e-6);
}
