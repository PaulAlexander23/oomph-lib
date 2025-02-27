#ifndef PARAMETER_FUNCTIONS_HEADER
#define PARAMETER_FUNCTIONS_HEADER

#include <functional>

#include "generic.h"
#include "generic/Vector.h"

namespace oomph
{
  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  slip_function_factory(const double& slip_length);

  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  wall_velocity_function_factory(double*& wall_velocity_pt);
} // namespace oomph
#endif
