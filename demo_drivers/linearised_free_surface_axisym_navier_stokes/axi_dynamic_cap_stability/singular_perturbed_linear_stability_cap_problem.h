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
  private:
    // List of augmented element numbers
    Vector<unsigned> Augmented_bulk_element_number;

  public:
    // Use the boundary id's from the base class
    // using PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
    //                                         PERTURBED_ELEMENT,
    //                                         TIMESTEPPER>::Boundary_id;
    enum Boundary_id
    {
      Upper_boundary_id,
      Outer_boundary_with_slip_id,
      Free_surface_boundary_id,
      Inner_boundary_id,
    };

    // Constructor
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

      augment_bulk_elements();

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

    void augment_bulk_elements()
    {
      double inner_radius = this->parameters_pt()->augmented_radius;

      // Ensure the two elements closest to the corner are augmented and get
      // their sizes
      PERTURBED_ELEMENT* corner_el_pt = 0;
      unsigned node_index = 0;

      double corner_element_size = 0.0;
      for (unsigned i = 0; i < 2; i++)
      {
        unsigned element_index;
        switch (i)
        {
          case 0:
            this->find_corner_bulk_element_and_node(
              Boundary_id::Outer_boundary_with_slip_id,
              Boundary_id::Free_surface_boundary_id,
              element_index,
              node_index);
            corner_el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
              this->fluid_mesh_pt()->boundary_element_pt(
                Boundary_id::Outer_boundary_with_slip_id, element_index));
            break;
          case 1:
            this->find_corner_bulk_element_and_node(
              Boundary_id::Free_surface_boundary_id,
              Boundary_id::Outer_boundary_with_slip_id,
              element_index,
              node_index);
            corner_el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
              this->fluid_mesh_pt()->boundary_element_pt(
                Boundary_id::Free_surface_boundary_id, element_index));
            break;
          default:
            break;
        }

        // corner_el_pt->augment();
        // corner_el_pt->add_additional_terms();
        // corner_el_pt->swap_unknowns();

        // Find the element iterator with the bulk mesh
        std::vector<GeneralisedElement*>::iterator iter =
          std::find(this->fluid_mesh_pt()->element_pt().begin(),
                    this->fluid_mesh_pt()->element_pt().end(),
                    dynamic_cast<GeneralisedElement*>(corner_el_pt));

        // Use this to get the element number
        unsigned e =
          std::distance(this->fluid_mesh_pt()->element_pt().begin(), iter);
        // Add the element number to the augmented element number vector
        Augmented_bulk_element_number.push_back(e);

        corner_element_size =
          std::max(corner_element_size, corner_el_pt->size());
      }

      if (inner_radius < 0)
      {
        inner_radius = 5.0 * pow(2.0 * corner_element_size, 0.5);
      }
      Node* contact_line_node_pt = corner_el_pt->node_pt(node_index);

      // Loop over the elements to set the consitutive law and jacobian
      unsigned n_bulk = this->fluid_mesh_pt()->nelement();
      for (unsigned e = 0; e < n_bulk; e++)
      {
        // Upcast from GeneralisedElement to the present element
        PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          this->fluid_mesh_pt()->element_pt(e));

        // Augmented elements close to the corner
        // Check distance from
        // s centre is centre of mass of a uniform triangle, so (1/3,1/3)
        // for Triangle[(1,0),(0,1),(0,0)]
        Vector<double> s_centre(2, 1.0 / 3.0);
        Vector<double> element_centre_x(2, 0.0);
        el_pt->get_x(s_centre, element_centre_x);
        double dist = 0;
        for (unsigned i = 0; i < 2; i++)
        {
          dist +=
            pow(element_centre_x[i] - contact_line_node_pt->position(i), 2.0);
        }
        dist = pow(dist, 0.5);

        // If the distance to the corner is within the "inner" region, ...
        if (dist < inner_radius)
        {
          // If this element is not already augmented, augment it
          // if (!el_pt->is_augmented())
          {
            // el_pt->augment();
            // el_pt->add_additional_terms();
            // el_pt->swap_unknowns();

            Augmented_bulk_element_number.push_back(e);
          }
        }
      }

      // Output the number of augmented elements
      oomph_info << Augmented_bulk_element_number.size()
                 << " augmented elements" << std::endl;
      if (Augmented_bulk_element_number.size() == 0)
      {
        oomph_info << "WARNING: No augmented elements! Try setting the "
                      "augmented region to be larger."
                   << std::endl;
      }
    }
  };
}; // namespace oomph
#endif
