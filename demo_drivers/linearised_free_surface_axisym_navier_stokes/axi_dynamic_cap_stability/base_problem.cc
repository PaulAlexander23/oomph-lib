#include "base_problem.h"

namespace oomph
{
  template class SingularAxisymDynamicCapProblem<
    SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>;
};
