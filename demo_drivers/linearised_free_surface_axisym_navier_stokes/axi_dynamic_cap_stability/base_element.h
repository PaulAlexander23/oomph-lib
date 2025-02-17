#ifndef BASE_ELEMENT_HEADER
#define BASE_ELEMENT_HEADER

#include "axisym_navier_stokes/singular_axisym_navier_stokes_elements.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/projectable_axisymmetric_Ttaylor_hood_elements.h"

namespace oomph
{
  extern template class SolidSingularAxisymNavierStokesElement<
    ProjectableAxisymmetricTTaylorHoodPVDElement>;
};

#endif
