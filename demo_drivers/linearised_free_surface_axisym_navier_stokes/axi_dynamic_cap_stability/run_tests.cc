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
  BOOST_TEST(abs(eigenvalue[0].real() - (-1.1117228111526063)) < 1e-6);
}


BOOST_AUTO_TEST_CASE(total_velocity_equations)
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

  perturbed_problem.set_always_take_one_newton_step();
  perturbed_problem.disable_singular_correction();
  // Pinning the mesh deformation
  perturbed_problem.pin_horizontal_mesh_deformation();
  perturbed_problem.pin_vertical_mesh_deformation();
  // Pin the momentum equations
  perturbed_problem.pin_volume_constraint();
  perturbed_problem.pin_fluid();
  perturbed_problem.set_constant_lagrange_free_surface_boundary_condition(0.0,
                                                                          0.0);
  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.pin_wall_velocity(upper, Vector<double>(6, 0.0));

  ofstream file("dofs.dat");
  perturbed_problem.describe_dofs(file);
  file.close();
  for (unsigned i = 0; i < 2; i++)
  {
    if (i == 0)
    {
      perturbed_problem.set_the_singular_correction(Vector<double>(2, 0.0));
    }
    else if (i == 1)
    {
      perturbed_problem.set_the_singular_correction(Vector<double>(2, 0.01));
    }
    else
    {
      perturbed_problem
        .unset_constant_lagrange_free_surface_boundary_condition();
    }


    perturbed_problem.setup_new_data();

    DoubleVector residuals;
    perturbed_problem.get_residuals(residuals);
    BOOST_TEST(abs(residuals.max()) < 1e-8);
    if (abs(residuals.max()) >= 1e-8)
    {
      residuals.output("residuals.dat");
    }

    // Steady newton solve
    try
    {
      perturbed_problem.steady_newton_solve();
      perturbed_problem.doc_solution();
    }
    catch (NewtonSolverError& error)
    {
      DoubleVector dummy_residuals;
      CRDoubleMatrix jacobian;
      perturbed_problem.get_jacobian(dummy_residuals, jacobian);
      jacobian.sparse_indexed_output("jacobian.dat", true);
    }
  }
}

BOOST_AUTO_TEST_CASE(internal_boundary_test_mode_0)
{
  Params parameters;
  parameters.slip_length = 0.0;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  base_problem.steady_newton_solve();
  base_problem.create_restart_file();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();
  // Pinning the mesh deformation
  perturbed_problem.pin_horizontal_mesh_deformation();
  perturbed_problem.pin_vertical_mesh_deformation();

  perturbed_problem.pin_volume_constraint();
  perturbed_problem.set_constant_lagrange_free_surface_boundary_condition(0.0,
                                                                          0.0);
  perturbed_problem.pin_flux_constraint();
  perturbed_problem.pin_wall_velocity(outer, Vector<double>(6, 0.0));
  perturbed_problem.pin_wall_velocity(upper, Vector<double>(6, 0.0));
  perturbed_problem.pin_wall_velocity(lower, Vector<double>(6, 0.0));
  perturbed_problem.pin_wall_velocity(inner, Vector<double>(6, 0.0));
  // The kinematic condition has lagrange multiplier contributions to the
  // momentum equations so should be imposed as a no penetration condition.
  perturbed_problem.set_always_take_one_newton_step();
  // Check the boundary conditions
  Vector<Vector<double>> wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  // Steady newton solve
  DoubleVector residuals;
  // The residuals should be zero given that there is no driving force away from
  // the equalibrium.
  perturbed_problem.get_residuals(residuals);
  BOOST_TEST(abs(residuals.max()) < 1e-8);


  DoubleVector dummy_residuals;
  CRDoubleMatrix jacobian;
  perturbed_problem.get_jacobian(dummy_residuals, jacobian);
  jacobian.sparse_indexed_output("jacobian.dat", true);
  ofstream file("dofs.dat");
  perturbed_problem.describe_dofs(file);
  file.close();

  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Set the singular solution contribution to something non-zero
  perturbed_problem.set_the_singular_correction(Vector<double>(2, 0.01));
  // And make sure the wall velocity is pinned to zero correctly, so the
  // new dirichlet conditions are used
  perturbed_problem.set_outer_boundary_condition();
  // Setup the new data for the augmented problem
  perturbed_problem.setup_new_data();

  // Check the wall dofs
  wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  perturbed_problem.doc_solution();
  DoubleVector new_residuals;
  // The residuals should be still be zero (or small) given that there is no
  // driving force away from the equalibrium.
  perturbed_problem.get_residuals(new_residuals);
  BOOST_TEST(abs(new_residuals.max()) < 1e-3);
  // The problem should also solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Check that the velocity is close to zero.
  Vector<Vector<double>> velocity = perturbed_problem.velocity();
  for (unsigned n = 0; n < velocity.size(); n++)
  {
    BOOST_TEST(abs(velocity[n][0]) < 1e-2);
    BOOST_TEST(abs(velocity[n][2]) < 1e-2);
  }
}

BOOST_AUTO_TEST_CASE(free_surface_fixed_c)
{
  Params parameters;
  parameters.azimuthal_mode_number = 1;
  parameters.slip_length = 0.0;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  base_problem.steady_newton_solve();
  base_problem.create_restart_file();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();
  // Pinning the mesh deformation
  // perturbed_problem.pin_horizontal_mesh_deformation();
  // perturbed_problem.pin_vertical_mesh_deformation();

  perturbed_problem.pin_volume_constraint();
  perturbed_problem.pin_wall_velocity(outer, Vector<double>(6, 0.0));
  perturbed_problem.pin_wall_velocity(upper, Vector<double>(6, 0.0));
  perturbed_problem.pin_centre_corner_lagrange_multipler();
  perturbed_problem.pin_contact_line_and_lagrange_multiplier();
  perturbed_problem.pin_wall_velocity(inner, Vector<double>(6, 0.0));
  // The kinematic condition has lagrange multiplier contributions to the
  // momentum equations so should be imposed as a no penetration condition.
  perturbed_problem.set_always_take_one_newton_step();
  // Check the boundary conditions
  Vector<Vector<double>> wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  // Steady newton solve
  DoubleVector residuals;
  // The residuals should be zero given that there is no driving force away from
  // the equalibrium.
  perturbed_problem.get_residuals(residuals);
  BOOST_TEST(abs(residuals.max()) < 1e-8);


  DoubleVector dummy_residuals;
  CRDoubleMatrix jacobian;
  perturbed_problem.get_jacobian(dummy_residuals, jacobian);
  jacobian.sparse_indexed_output("jacobian.dat", true);
  ofstream file("dofs.dat");
  perturbed_problem.describe_dofs(file);
  file.close();

  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Set the singular solution contribution to something non-zero
  perturbed_problem.set_the_singular_correction(Vector<double>(2, 0.01));
  // And make sure the wall velocity is pinned to zero correctly, so the
  // new dirichlet conditions are used
  perturbed_problem.set_outer_boundary_condition();
  // Setup the new data for the augmented problem
  perturbed_problem.setup_new_data();

  // Check the wall dofs
  wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  perturbed_problem.doc_solution();
  DoubleVector new_residuals;
  // The residuals should be still be zero (or small) given that there is no
  // driving force away from the equalibrium.
  perturbed_problem.get_residuals(new_residuals);
  BOOST_TEST(abs(new_residuals.max()) < 1e-2);
  // The problem should also solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Check that the velocity is close to zero.
  Vector<Vector<double>> velocity = perturbed_problem.velocity();
  for (unsigned n = 0; n < velocity.size(); n++)
  {
    BOOST_TEST(abs(velocity[n][0]) < 1e-4);
    BOOST_TEST(abs(velocity[n][2]) < 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(dirichlet_bcs)
{
  Params parameters;
  parameters.azimuthal_mode_number = 0;
  parameters.slip_length = 0.0;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  base_problem.steady_newton_solve();
  base_problem.create_restart_file();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();
  // Pinning the mesh deformation
  // perturbed_problem.pin_horizontal_mesh_deformation();
  // perturbed_problem.pin_vertical_mesh_deformation();

  perturbed_problem.pin_volume_constraint();
  perturbed_problem.pin_wall_velocity(upper, Vector<double>(6, 0.0));
  perturbed_problem.pin_centre_corner_lagrange_multipler();
  perturbed_problem.pin_contact_line_and_lagrange_multiplier();
  perturbed_problem.pin_wall_velocity(inner, Vector<double>(6, 0.0));
  // The kinematic condition has lagrange multiplier contributions to the
  // momentum equations so should be imposed as a no penetration condition.
  perturbed_problem.set_always_take_one_newton_step();
  // Check the boundary conditions
  Vector<Vector<double>> wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  // Steady newton solve
  DoubleVector residuals;
  // The residuals should be zero given that there is no driving force away from
  // the equalibrium.
  perturbed_problem.get_residuals(residuals);
  BOOST_TEST(abs(residuals.max()) < 1e-8);


  DoubleVector dummy_residuals;
  CRDoubleMatrix jacobian;
  perturbed_problem.get_jacobian(dummy_residuals, jacobian);
  jacobian.sparse_indexed_output("jacobian.dat", true);
  ofstream file("dofs.dat");
  perturbed_problem.describe_dofs(file);
  file.close();

  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Set the singular solution contribution to something non-zero
  perturbed_problem.set_the_singular_correction(Vector<double>(2, 0.01));
  // And make sure the wall velocity is pinned to zero correctly, so the
  // new dirichlet conditions are used
  perturbed_problem.set_outer_boundary_condition();
  // Setup the new data for the augmented problem
  perturbed_problem.setup_new_data();

  // Check the wall dofs
  wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  perturbed_problem.doc_solution();
  DoubleVector new_residuals;
  // The residuals should be still be zero (or small) given that there is no
  // driving force away from the equalibrium.
  perturbed_problem.get_residuals(new_residuals);
  BOOST_TEST(abs(new_residuals.max()) < 1e-2);
  // The problem should also solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Check that the velocity is close to zero.
  Vector<Vector<double>> velocity = perturbed_problem.velocity();
  for (unsigned n = 0; n < velocity.size(); n++)
  {
    BOOST_TEST(abs(velocity[n][0]) < 1e-4);
    BOOST_TEST(abs(velocity[n][2]) < 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(smooth_velocity_on_outer_wall)
{
  Params parameters;
  parameters.azimuthal_mode_number = 0;
  parameters.slip_length = 0.0;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  base_problem.steady_newton_solve();
  base_problem.create_restart_file();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();
  // Pinning the mesh deformation
  perturbed_problem.pin_horizontal_mesh_deformation();
  perturbed_problem.pin_vertical_mesh_deformation();

  perturbed_problem.pin_wall_velocity(lower, Vector<double>(6, 0.0));
  perturbed_problem.pin_wall_velocity(upper, Vector<double>(6, 0.0));

  perturbed_problem.set_constant_lagrange_free_surface_boundary_condition(0.0,
                                                                          0.0);

  perturbed_problem.pin_volume_constraint();
  perturbed_problem.pin_flux_constraint();

  perturbed_problem.pin_wall_velocity(inner, Vector<double>(6, 0.0));

  // The kinematic condition has lagrange multiplier contributions to the
  // momentum equations so should be imposed as a no penetration condition.
  perturbed_problem.set_always_take_one_newton_step();

  // Check the boundary conditions
  Vector<Vector<double>> wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  // Steady newton solve
  DoubleVector residuals;
  // The residuals should be zero given that there is no driving force away from
  // the equalibrium.
  perturbed_problem.get_residuals(residuals);
  BOOST_TEST(abs(residuals.max()) < 1e-8);


  DoubleVector dummy_residuals;
  CRDoubleMatrix jacobian;
  perturbed_problem.get_jacobian(dummy_residuals, jacobian);
  jacobian.sparse_indexed_output("jacobian.dat", true);
  ofstream file("dofs.dat");
  perturbed_problem.describe_dofs(file);
  file.close();

  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Set the singular solution contribution to something non-zero
  perturbed_problem.set_the_singular_correction(Vector<double>(2, 0.01));
  // And make sure the wall velocity is pinned to zero correctly, so the
  // new dirichlet conditions are used
  perturbed_problem.set_outer_boundary_condition();
  // Setup the new data for the augmented problem
  perturbed_problem.setup_new_data();

  // Check the wall dofs
  wall_velocity = perturbed_problem.wall_velocity();
  for (unsigned n = 0; n < wall_velocity.size(); n++)
  {
    BOOST_TEST(abs(wall_velocity[n][1] - 0.0) < 1e-8);
  }

  perturbed_problem.doc_solution();
  DoubleVector new_residuals;
  // The residuals should be still be zero (or small) given that there is no
  // driving force away from the equalibrium.
  perturbed_problem.get_residuals(new_residuals);
  BOOST_TEST(abs(new_residuals.max()) < 1e-3);
  // The problem should also solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Check that the velocity is close to zero.
  Vector<Vector<double>> velocity = perturbed_problem.velocity();
  for (unsigned n = 0; n < velocity.size(); n++)
  {
    BOOST_TEST(abs(velocity[n][0]) < 1e-2);
    BOOST_TEST(abs(velocity[n][2]) < 1e-2);
  }
}

// Fix the singular scaling to test the internal and external boundary
// conditions
BOOST_AUTO_TEST_CASE(mode_0_fix_c)
{
  Params parameters;
  // parameters.restart_filename = "RESLT/restart0.dat";
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  //  ifstream restart_filestream;
  //  restart_filestream.open(parameters.restart_filename);
  //  bool is_unsteady_restart = false;
  //  base_problem.read(restart_filestream, is_unsteady_restart);
  //  restart_filestream.close();
  base_problem.steady_newton_solve();
  base_problem.create_restart_file();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.disable_singular_correction();

  // Steady newton solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Eigensolve
  Vector<std::complex<double>> eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);

  perturbed_problem.set_the_singular_correction(Vector<double>(2, 1));
  perturbed_problem.setup_new_data();
  perturbed_problem.set_boundary_conditions();

  DoubleVector residuals;
  perturbed_problem.get_residuals(residuals);
  residuals.output("residuals.dat");
  ofstream file("dofs.dat");
  perturbed_problem.describe_dofs(file);
  file.close();
  perturbed_problem.doc_solution();

  // Eigensolve
  eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);

  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Eigensolve
  eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
}

// Fix the singular scaling to test the internal and external boundary
// conditions
// BOOST_AUTO_TEST_CASE(mode_0)
// {
//   Params parameters;
//   // parameters.restart_filename = "RESLT/restart0.dat";
//   parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;
// 
//   // Create the base problem
//   BASE_PROBLEM base_problem(&parameters);
//   //  ifstream restart_filestream;
//   //  restart_filestream.open(parameters.restart_filename);
//   //  bool is_unsteady_restart = false;
//   //  base_problem.read(restart_filestream, is_unsteady_restart);
//   //  restart_filestream.close();
//   base_problem.steady_newton_solve();
//   base_problem.create_restart_file();
// 
//   base_problem.reset_lagrange();
//   base_problem.assign_initial_values_impulsive();
// 
//   // Create the linear problem
//   PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
//                                       base_problem.free_surface_mesh_pt(),
//                                       base_problem.slip_surface_mesh_pt(),
//                                       &parameters);
// 
//   perturbed_problem.assign_initial_values_impulsive();
//   perturbed_problem.set_always_take_one_newton_step();
// 
//   DoubleVector dummy_residuals;
//   CRDoubleMatrix jacobian;
//   perturbed_problem.get_jacobian(dummy_residuals, jacobian);
//   jacobian.sparse_indexed_output("jacobian.dat", true);
// 
//   ofstream file("dofs.dat");
//   perturbed_problem.describe_dofs(file);
//   file.close();
// 
//   // debug_jacobian(&perturbed_problem);
// 
//   // Steady newton solve
//   perturbed_problem.steady_newton_solve();
//   perturbed_problem.doc_solution();
// 
//   // Eigensolve
//   Vector<std::complex<double>> eigenvalue =
//     perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
// }

// Fix the singular scaling to test the internal and external boundary
// conditions
//BOOST_AUTO_TEST_CASE(mode_1_fix_c)
//{
//  Params parameters;
//  // parameters.restart_filename = "RESLT/restart0.dat";
//  parameters.azimuthal_mode_number = 1;
//  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;
//
//  // Create the base problem
//  BASE_PROBLEM base_problem(&parameters);
//  //  ifstream restart_filestream;
//  //  restart_filestream.open(parameters.restart_filename);
//  //  bool is_unsteady_restart = false;
//  //  base_problem.read(restart_filestream, is_unsteady_restart);
//  //  restart_filestream.close();
//  base_problem.steady_newton_solve();
//  base_problem.create_restart_file();
//
//  base_problem.reset_lagrange();
//  base_problem.assign_initial_values_impulsive();
//
//  // Create the linear problem
//  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
//                                      base_problem.free_surface_mesh_pt(),
//                                      base_problem.slip_surface_mesh_pt(),
//                                      &parameters);
//
//  perturbed_problem.assign_initial_values_impulsive();
//  perturbed_problem.disable_singular_correction();
//
//  // Steady newton solve
//  perturbed_problem.steady_newton_solve();
//  perturbed_problem.doc_solution();
//
//  // Eigensolve
//  Vector<std::complex<double>> eigenvalue =
//    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
//
//  perturbed_problem.set_the_singular_correction(Vector<double>(2, 0.1));
//  perturbed_problem.setup_new_data();
//  perturbed_problem.set_boundary_conditions();
//
//  DoubleVector residuals;
//  perturbed_problem.get_residuals(residuals);
//  residuals.output("residuals.dat");
//  ofstream file("dofs.dat");
//  perturbed_problem.describe_dofs(file);
//  file.close();
//  perturbed_problem.doc_solution();
//
//  // Eigensolve
//  eigenvalue =
//    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
//
//  perturbed_problem.steady_newton_solve();
//  perturbed_problem.doc_solution();
//
//  // Eigensolve
//  eigenvalue =
//    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
//}

// Fix the singular scaling to test the internal and external boundary
// conditions
BOOST_AUTO_TEST_CASE(mode_1)
{
  Params parameters;
  // parameters.restart_filename = "RESLT/restart0.dat";
  parameters.azimuthal_mode_number = 1;
  parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;

  // Create the base problem
  BASE_PROBLEM base_problem(&parameters);
  //  ifstream restart_filestream;
  //  restart_filestream.open(parameters.restart_filename);
  //  bool is_unsteady_restart = false;
  //  base_problem.read(restart_filestream, is_unsteady_restart);
  //  restart_filestream.close();
  base_problem.steady_newton_solve();
  base_problem.create_restart_file();

  base_problem.reset_lagrange();
  base_problem.assign_initial_values_impulsive();

  // Create the linear problem
  PERTURBED_PROBLEM perturbed_problem(base_problem.bulk_mesh_pt(),
                                      base_problem.free_surface_mesh_pt(),
                                      base_problem.slip_surface_mesh_pt(),
                                      &parameters);

  perturbed_problem.assign_initial_values_impulsive();
  perturbed_problem.set_always_take_one_newton_step();

  DoubleVector dummy_residuals;
  CRDoubleMatrix jacobian;
  perturbed_problem.get_jacobian(dummy_residuals, jacobian);
  jacobian.sparse_indexed_output("jacobian.dat", true);

  ofstream file("dofs.dat");
  perturbed_problem.describe_dofs(file);
  file.close();

  // debug_jacobian(&perturbed_problem);

  // Steady newton solve
  perturbed_problem.steady_newton_solve();
  perturbed_problem.doc_solution();

  // Eigensolve
  Vector<std::complex<double>> eigenvalue =
    perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
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

//// Augmented linear problem
// BOOST_AUTO_TEST_CASE(augmented_linear_problem)
//{
//   // Create the parameters
//   Params parameters;
//   parameters.azimuthal_mode_number = 1;
//   parameters.contact_angle = 120.0 / 180.0 * MathematicalConstants::Pi;
//   //*parameters.wall_velocity_pt = 0.01;
//
//   // Create the base problem
//   typedef SolidSingularAxisymNavierStokesElement<
//     ProjectableAxisymmetricTTaylorHoodPVDElement>
//     BASE_ELEMENT;
//   typedef BDF<2> TIMESTEPPER;
//   SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER> base_problem(
//     &parameters);
//   base_problem.steady_newton_solve();
//   base_problem.reset_lagrange();
//   base_problem.assign_initial_values_impulsive();
//
//   parameters.azimuthal_mode_number = 0;
//   unsigned local_doc_number = 0;
//
//   // Create the linear problem
//   typedef SingularOverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
//   SingularPerturbedLinearStabilityCapProblem<BASE_ELEMENT,
//                                              PERTURBED_ELEMENT,
//                                              TIMESTEPPER>
//     perturbed_problem(base_problem.bulk_mesh_pt(),
//                       base_problem.free_surface_mesh_pt(),
//                       base_problem.slip_surface_mesh_pt(),
//                       &parameters);
//
//   perturbed_problem.assign_initial_values_impulsive();
//   perturbed_problem.doc_solution();
//
//   // Steady newton solve
//   perturbed_problem.steady_newton_solve();
//   perturbed_problem.doc_solution();
//
//   // perturbed_problem.perturb_vertical_velocity();
//   // perturbed_problem.assign_initial_values_impulsive();
//   // perturbed_problem.doc_solution();
//
//   //// Time step (docs within the timestepper)
//   // perturbed_problem.timestep(0.01, 0.01);
//
//   // perturbed_problem.make_unsteady();
//   perturbed_problem.pin_horizontal_mesh_deformation();
//
//
//   // Eigensolve
//   Vector<std::complex<double>> eigenvalue =
//     perturbed_problem.solve_and_document_n_most_unstable_eigensolutions(1);
//   BOOST_TEST(abs(eigenvalue[0].real() - (-1.5677282757561151)) < 1e-6);
//
//   local_doc_number = perturbed_problem.doc_info().number();
//
//   // debug_jacobian(&perturbed_problem);
//
//   parameters.azimuthal_mode_number = 1;
//
//   // Create the linear problem
//   SingularPerturbedLinearStabilityCapProblem<BASE_ELEMENT,
//                                              PERTURBED_ELEMENT,
//                                              TIMESTEPPER>
//     perturbed_problem0(base_problem.bulk_mesh_pt(),
//                        base_problem.free_surface_mesh_pt(),
//                        base_problem.slip_surface_mesh_pt(),
//                        &parameters);
//   perturbed_problem0.assign_initial_values_impulsive();
//   perturbed_problem0.doc_info().number() = local_doc_number;
//   perturbed_problem0.steady_newton_solve();
//   perturbed_problem0.make_unsteady();
//   perturbed_problem0.pin_horizontal_mesh_deformation();
//   eigenvalue =
//     perturbed_problem0.solve_and_document_n_most_unstable_eigensolutions(1);
//   BOOST_TEST(abs(eigenvalue[0].real() - (-0.59043964949883476)) < 1e-6);
// }
