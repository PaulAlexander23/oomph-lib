#include <boost/test/unit_test.hpp>

#include "generic.h"
#include "utility_functions.h"

using namespace oomph;

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

BOOST_AUTO_TEST_CASE(compare_matrix_different)
{
  const unsigned N = 4;
  DenseDoubleMatrix A(N, N, 0.0);
  DenseDoubleMatrix B(N, N, 0.0);
  A(2, 1) = 1.0;
  B(2, 1) = 1e-10;
  BOOST_TEST(compare_matrices(A, B) == 0);
}
