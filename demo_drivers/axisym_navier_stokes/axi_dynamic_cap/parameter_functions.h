#ifndef PARAMETER_FUNCTIONS_HEADER
#define PARAMETER_FUNCTIONS_HEADER

#include "generic.h"

namespace oomph
{
  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  slip_function_factory(const double& slip_length)
  {
    // Return a lambda that captures 's' and takes 't' and 'x' as arguments
    return [slip_length](const double& t,
                         const Vector<double>& x,
                         const Vector<double>& n,
                         Vector<double>& slip) -> void
    {
      slip[0] = -1;
      slip[1] = slip_length;
      slip[2] = -1;
    };
  }

  std::function<void(const double&,
                     const Vector<double>&,
                     const Vector<double>&,
                     Vector<double>&)>
  wall_velocity_function_factory(double*& wall_velocity_pt)
  {
    // Return a lambda that captures 's' and takes 't' and 'x' as arguments
    return [&wall_velocity_pt](const double& t,
                               const Vector<double>& x,
                               const Vector<double>& n,
                               Vector<double>& velocity) -> void
    {
      // Assign solution
      velocity[0] = 0.0;
      velocity[1] = *wall_velocity_pt;
      velocity[2] = 0.0;
    };
  }
} // namespace oomph
#endif
