#ifndef PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER
#define PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER

#include "perturbed_linear_stability_cap_problem_base.h"

namespace oomph
{
  //===========start_of_perturbed_linear_stability_cap_problem_class============
  // A class that solves the linearised Axisymmetric Navier--Stokes equations
  // to compute the stability of a dynamic interface in a cylindrical container
  // with imposed contact angle at the boundary.
  //============================================================================
  template<class BASE_ELEMENT, class PERTURBED_ELEMENT, class TIMESTEPPER>
  class PerturbedLinearStabilityCapProblem
    : public PerturbedLinearStabilityCapProblemBase<BASE_ELEMENT,
                                                    PERTURBED_ELEMENT,
                                                    TIMESTEPPER>
  {
  public:
    // Constructor
    // Uses the external base and surface mesh of the base state to create a
    // new problem solving the linearised system on the same mesh geometry.
    // This was easier than using the multi domain elements due to the
    // complexities on the free surface boundary and the contact line.
    // The azimuthal_mode_number is passed in as a parameter.
    PerturbedLinearStabilityCapProblem(Mesh* external_base_mesh_pt,
                                       Mesh* external_free_surface_mesh_pt,
                                       Mesh* external_slip_surface_mesh_pt,
                                       Params* const& params_pt)
      : PerturbedLinearStabilityCapProblemBase<BASE_ELEMENT,
                                               PERTURBED_ELEMENT,
                                               TIMESTEPPER>(
          external_base_mesh_pt,
          external_free_surface_mesh_pt,
          external_slip_surface_mesh_pt,
          params_pt)
    {
    }
  };
} // namespace oomph

#endif
