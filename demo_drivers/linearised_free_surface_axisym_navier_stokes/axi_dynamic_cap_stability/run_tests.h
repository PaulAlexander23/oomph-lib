#ifndef RUN_TESTS_HEADER
#define RUN_TESTS_HEADER

// STD includes
#include <iostream>

// OOMPH-LIB include files
#include "generic/timesteppers.h"
// Other demo includes
#include "../../axisym_navier_stokes/axi_dynamic_cap/projectable_axisymmetric_Ttaylor_hood_elements.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/singular_axisym_dynamic_cap_problem.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/parameters.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/utility_functions.h"
// Local includes
#include "base_element.h"
#include "singular_perturbed_linear_stability_cap_problem.h"
#include "singular_overlaying_my_linear_elements.h"

namespace oomph
{
  typedef SolidSingularAxisymNavierStokesElement<
    ProjectableAxisymmetricTTaylorHoodPVDElement>
    BASE_ELEMENT;
  typedef SingularOverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  typedef SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>
    BASE_PROBLEM;
  typedef SingularPerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                                     PERTURBED_ELEMENT,
                                                     TIMESTEPPER>
    PERTURBED_PROBLEM;

  enum
  {
    upper,
    outer,
    lower,
    inner,
  };

} // namespace oomph

#endif
