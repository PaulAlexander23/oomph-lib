#ifndef SINGULAR_PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER
#define SINGULAR_PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER

#include "perturbed_linear_stability_cap_problem.h"

namespace oomph
{
  template<class BASE_ELEMENT, class PERTURBED_ELEMENT, class TIMESTEPPER>
  class SingularPerturbedLinearStabilityCapProblem
    : public PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                                PERTURBED_ELEMENT,
                                                TIMESTEPPER>
  {
  public:
    SingularPerturbedLinearStabilityCapProblem(
      Mesh* external_base_mesh_pt,
      Mesh* external_free_surface_mesh_pt,
      Mesh* external_slip_surface_mesh_pt,
      Params* const& params_pt)
      : PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                           PERTURBED_ELEMENT,
                                           TIMESTEPPER>(
          external_base_mesh_pt,
          external_free_surface_mesh_pt,
          external_slip_surface_mesh_pt,
          params_pt)
    {
      PerturbedLinearStabilityCapProblem<
        BASE_ELEMENT,
        PERTURBED_ELEMENT,
        TIMESTEPPER>::remove_boundary_elements();

      this->add_boundary_elements();

      // Set up the connections to the base state
      this->set_up_overlapping_domain_functions();

      // Set the boundary conditions
      this->set_boundary_conditions();

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Set up the equation numbering so we are ready to solve the problem.
      oomph_info << "Number of unknowns: " << this->assign_eqn_numbers()
                 << std::endl;
    }
  };
}; // namespace oomph
#endif
