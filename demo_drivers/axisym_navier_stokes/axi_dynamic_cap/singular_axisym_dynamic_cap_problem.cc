#include "singular_axisym_dynamic_cap_problem.h"

namespace oomph
{
  template class SingularAxisymDynamicCapProblem<
    SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>;
};
