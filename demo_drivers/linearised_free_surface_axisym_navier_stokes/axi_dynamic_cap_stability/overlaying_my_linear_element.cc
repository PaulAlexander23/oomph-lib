#include "overlaying_my_linear_element.h"

namespace oomph
{
  template class OverlayingMyLinearElement<
    SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>>;
  template class FaceGeometry<
    OverlayingMyLinearElement<SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>>>;
}; // namespace oomph
