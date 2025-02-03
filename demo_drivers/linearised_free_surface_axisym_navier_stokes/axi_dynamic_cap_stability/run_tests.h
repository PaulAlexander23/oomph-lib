#ifndef RUN_TESTS_HEADER
#define RUN_TESTS_HEADER

// STD includes
#include <iostream>

// OOMPH-LIB include files
#include "generic.h"
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "constitutive.h"
#include "solid.h"
#include "meshes/triangle_mesh.h"
#include "linearised_axisym_navier_stokes.h"

// Other demo includes
#include "../../axisym_navier_stokes/axi_dynamic_cap/projectable_axisymmetric_Ttaylor_hood_elements.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/singular_axisym_dynamic_cap_problem.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/parameters.h"
#include "../../axisym_navier_stokes/axi_dynamic_cap/utility_functions.h"

// Local includes
#include "linearised_axisymmetric_fluid_interface_elements.h"
#include "decomposed_linear_elasticity_elements.h"
#include "linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_linearised_elastic_axisym_fluid_interface_element.h"
#include "overlaying_my_linear_element.h"
#include "perturbed_linear_stability_cap_problem.h"

namespace oomph
{
  typedef SingularAxisymNavierStokesElement<
    ProjectableAxisymmetricTTaylorHoodPVDElement>
    BASE_ELEMENT;
  // typedef OverlayingTLinearisedAxisymNSPVDElement PERTURBED_ELEMENT;
  typedef OverlayingMyLinearElement<BASE_ELEMENT> PERTURBED_ELEMENT;
  typedef BDF<2> TIMESTEPPER;
  typedef SingularAxisymDynamicCapProblem<BASE_ELEMENT, TIMESTEPPER>
    AXISYM_PROBLEM;
  typedef PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                             PERTURBED_ELEMENT,
                                             TIMESTEPPER>
    PERTURBED_PROBLEM;

  std::shared_ptr<AXISYM_PROBLEM> createBaseProblem();
  std::shared_ptr<PERTURBED_PROBLEM> createLinearProblem(
    std::shared_ptr<AXISYM_PROBLEM> base_problem_pt, Params& parameters);

  enum
  {
    upper,
    outer,
    lower,
    inner,
  };

} // namespace oomph

#endif
