#include "singular_perturbed_linear_stability_cap_problem.h"

namespace oomph
{
  // Explicit instantiation
  template class SingularPerturbedLinearStabilityCapProblem<
    SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    SingularOverlayingMyLinearElement<SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>>,
    BDF<2>>;
}; // namespace oomph
