#ifndef BASE_PROBLEM_HEADER
#define BASE_PROBLEM_HEADER

#include "generic/timesteppers.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/projectable_axisymmetric_Ttaylor_hood_elements.h"
#include "axisym_navier_stokes/singular_axisym_navier_stokes_elements.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/singular_axisym_dynamic_cap_problem.h"

namespace oomph
{
  extern template class SingularAxisymDynamicCapProblem<
    SolidSingularAxisymNavierStokesElement<
      ProjectableAxisymmetricTTaylorHoodPVDElement>,
    BDF<2>>;
};


#endif
