#define BOOST_TEST_MODULE grad_velocity_test_module
#include <boost/test/included/unit_test.hpp>

#include "eigensolution_functions.h"

using namespace std;
using namespace oomph;
using namespace oomph::parameters;

BOOST_AUTO_TEST_CASE(grad_velocity)
{
  Vector<double> x_centre(2, 0.0);
  x_centre[0] = 1.0;
  x_centre_pt = &x_centre;

  Vector<Vector<double>> expected_grad(2);
  for (unsigned d = 0; d < 2; d++)
  {
    expected_grad[d].resize(2);
    for (unsigned k = 0; k < 2; k++)
    {
      expected_grad[d][k] = 0.0;
    }
  }

  Vector<double> x0(2, 0.0);

  x0[0] = 0.9;
  x0[1] = 0.1;
  Vector<double> u0 = axisym_velocity_singular_fct(x0);
  const double dx = 1e-8;
  for (unsigned d = 0; d < 2; d++)
  {
    Vector<double> x1 = x0;
    x1[d] -= dx;
    Vector<double> u1 = axisym_velocity_singular_fct(x1);
    for (unsigned k = 0; k < 2; k++)
    {
      expected_grad[k][d] = (u0[k] - u1[k]) / dx;
    }
  }

  Vector<Vector<double>> actual_grad = grad_axisym_velocity_singular_fct(x0);

  for (unsigned d = 0; d < 2; d++)
  {
    for (unsigned k = 0; k < 2; k++)
    {
      BOOST_TEST(abs(actual_grad[d][k] - expected_grad[d][k]) < 4e-8);
    }
  }
}
