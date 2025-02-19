#include "singular_overlaying_my_linear_elements.h"
#include "base_element.h"

namespace oomph
{
  template class SingularOverlayingMyLinearElement<
    SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>>;

  template class FaceGeometry<SingularOverlayingMyLinearElement<SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>>>;
};
