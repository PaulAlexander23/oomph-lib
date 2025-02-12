#ifndef SINGULAR_PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER
#define SINGULAR_PERTURBED_LINEAR_STABILITY_CAP_PROBLEM_HEADER

#include "perturbed_linear_stability_cap_problem.h"
#include "decomposed_pressure_evaluation_elements.h"

namespace oomph
{
  template<class BASE_ELEMENT, class PERTURBED_ELEMENT, class TIMESTEPPER>
  class SingularPerturbedLinearStabilityCapProblem
    : public PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
                                                PERTURBED_ELEMENT,
                                                TIMESTEPPER>
  {
  public:
    typedef SingularNavierStokesSolutionElement<
      OverlayingMyLinearElement<BASE_ELEMENT>>
      SCALING_ELEMENT;

  private:
    // List of augmented element numbers
    Vector<unsigned> Augmented_bulk_element_number;

    // Singular solution scaling mesh
    Mesh* Singularity_scaling_mesh_pt;

    // Pressure contribution meshes
    Mesh* Pressure_contribution_mesh_1_pt;
    Mesh* Pressure_contribution_mesh_2_pt;

    // Eigensolution functions
    std::function<Vector<double>(const Vector<double>&)>
      Velocity_singular_function;
    std::function<Vector<Vector<double>>(const Vector<double>&)>
      Grad_velocity_singular_function;
    Node* Contact_line_node_pt;

  public:
    // Boundary ids enumeration
    enum Boundary_id
    {
      Upper_boundary_id,
      Outer_boundary_with_slip_id,
      Free_surface_boundary_id,
      Inner_boundary_id,
    };
    // Can't seem to use the boundary id's from the base class
    // using PerturbedLinearStabilityCapProblem<BASE_ELEMENT,
    //                                         PERTURBED_ELEMENT,
    //                                         TIMESTEPPER>::Boundary_id;

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
          params_pt),
        Singularity_scaling_mesh_pt(0),
        Pressure_contribution_mesh_1_pt(0),
        Pressure_contribution_mesh_2_pt(0)
    {
      // Setup the singular functions
      Contact_line_node_pt = this->find_corner_node(Outer_boundary_with_slip_id,
                                                    Free_surface_boundary_id);
      Velocity_singular_function = velocity_singular_function_factory(
        this->parameters_pt()->contact_angle, Contact_line_node_pt);
      Grad_velocity_singular_function = grad_velocity_singular_function_factory(
        this->parameters_pt()->contact_angle, Contact_line_node_pt);


      // Remove the original problem's boundary elements
      PerturbedLinearStabilityCapProblem<
        BASE_ELEMENT,
        PERTURBED_ELEMENT,
        TIMESTEPPER>::remove_boundary_elements();

      // Augment the bulk elements
      augment_bulk_elements();

      // Add the new sub meshes
      Singularity_scaling_mesh_pt = new Mesh;
      this->add_sub_mesh(Singularity_scaling_mesh_pt);
      Pressure_contribution_mesh_1_pt = new Mesh;
      this->add_sub_mesh(Pressure_contribution_mesh_1_pt);
      Pressure_contribution_mesh_2_pt = new Mesh;
      this->add_sub_mesh(Pressure_contribution_mesh_2_pt);

      this->add_boundary_elements();

      create_singularity_scaling_elements();
      create_pressure_contribution_1_elements();
      create_pressure_contribution_2_elements();

      // Set up the connections to the base state
      this->set_up_overlapping_domain_functions();

      setup_mesh_interaction();

      // Set the boundary conditions
      this->set_boundary_conditions();

      // Rebuild the global mesh
      this->rebuild_global_mesh();

      // Set up the equation numbering so we are ready to solve the problem.
      oomph_info << "Number of unknowns: " << this->assign_eqn_numbers()
                 << std::endl;
    }

    /// Augment the bulk elements within a small radius of the corner
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

        corner_el_pt->augment();

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
          if (!el_pt->is_augmented())
          {
            el_pt->augment();

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

    /// Create the singular solution scaling elements
    void create_singularity_scaling_elements()
    {
      oomph_info << "create_singularity_scaling_elements" << std::endl;
      // Create two scaling elements
      for (unsigned i = 0; i < 2; i++)
      {
        SCALING_ELEMENT* el_pt = new SCALING_ELEMENT;

        // Set the pointer to the velocity singular function for this
        // element, defined in parameters namespace
        el_pt->velocity_singular_fct() = Velocity_singular_function;

        // Set the pointer to the gradient of the velocity singular
        // function for this element, defined in parameters namespace
        el_pt->grad_velocity_singular_fct() = Grad_velocity_singular_function;

        // Set the pointer to the first pressure singular function for this
        // element, defined in parameters namespace
        el_pt->pressure_singular_fct_pt() = &pressure_singular_fct;

        // The singular function satisfies the Stokes equation
        el_pt->singular_function_satisfies_stokes_equation() = false;

        el_pt->pin_c();
        el_pt->set_c(0.0);

        // Add element to the mesh
        Singularity_scaling_mesh_pt->add_element_pt(el_pt);
      }
    }

    /// Create the pressure contribution elements for the first boundary
    void create_pressure_contribution_1_elements()
    {
      oomph_info << "create_pressure_contribution_1_elements" << std::endl;

      PERTURBED_ELEMENT* element_pt = 0;
      int face_index = 0;
      find_corner_bulk_element_and_face_index(Outer_boundary_with_slip_id,
                                              Free_surface_boundary_id,
                                              element_pt,
                                              face_index);


      const unsigned pressure_value_index = 0;
      DecomposedPressureEvaluationElement<PERTURBED_ELEMENT>* el_pt =
        new DecomposedPressureEvaluationElement<PERTURBED_ELEMENT>(
          element_pt, face_index, Contact_line_node_pt, pressure_value_index);

      // Add the two singularity solution scaling data
      el_pt->add_scaling_data(
        Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));
      el_pt->add_scaling_data(
        Singularity_scaling_mesh_pt->element_pt(1)->internal_data_pt(0));

      el_pt->set_boundary_number_in_bulk_mesh(Outer_boundary_with_slip_id);
      // Set the product of the Reynolds number and the inverse of the
      // Froude number
      el_pt->re_invfr_pt() =
        this->parameters_pt()->reynolds_inverse_froude_number_pt;
      // Set the direction of gravity
      el_pt->g_pt() = &this->parameters_pt()->gravity_vector;

      Pressure_contribution_mesh_1_pt->add_element_pt(el_pt);
    }

    /// Create the pressure contribution elements for the second boundary
    void create_pressure_contribution_2_elements()
    {
      oomph_info << "create_pressure_contribution_1_elements" << std::endl;

      PERTURBED_ELEMENT* element_pt = 0;
      int face_index = 0;
      find_corner_bulk_element_and_face_index(Free_surface_boundary_id,
                                              Outer_boundary_with_slip_id,
                                              element_pt,
                                              face_index);


      const unsigned pressure_value_index = 0;
      DecomposedPressureEvaluationElement<PERTURBED_ELEMENT>* el_pt =
        new DecomposedPressureEvaluationElement<PERTURBED_ELEMENT>(
          element_pt, face_index, Contact_line_node_pt, pressure_value_index);

      // Add the two singularity solution scaling data
      el_pt->add_scaling_data(
        Singularity_scaling_mesh_pt->element_pt(0)->internal_data_pt(0));
      el_pt->add_scaling_data(
        Singularity_scaling_mesh_pt->element_pt(1)->internal_data_pt(0));

      el_pt->set_boundary_number_in_bulk_mesh(Free_surface_boundary_id);
      // Set the product of the Reynolds number and the inverse of the
      // Froude number
      el_pt->re_invfr_pt() =
        this->parameters_pt()->reynolds_inverse_froude_number_pt;
      // Set the direction of gravity
      el_pt->g_pt() = &this->parameters_pt()->gravity_vector;
      el_pt->set_subtract_from_residuals();

      Pressure_contribution_mesh_2_pt->add_element_pt(el_pt);
    }

    /// Find the corner element and the face for the first of the two boundaries
    void find_corner_bulk_element_and_face_index(const unsigned& boundary_1_id,
                                                 const unsigned& boundary_2_id,
                                                 PERTURBED_ELEMENT*& element_pt,
                                                 int& face_index)
    {
      unsigned n_boundary_element =
        this->fluid_mesh_pt()->nboundary_element(boundary_1_id);
      for (unsigned e = 0; e < n_boundary_element; e++)
      {
        // Locally cache the element pointer
        FiniteElement* bulk_el_pt =
          this->fluid_mesh_pt()->boundary_element_pt(boundary_1_id, e);

        // Read out number of nodes in the element
        unsigned n_node = bulk_el_pt->nnode();
        for (unsigned i_node = 0; i_node < n_node; i_node++)
        {
          // If the node is on the free surface boundary as well then ...
          if (bulk_el_pt->node_pt(i_node)->is_on_boundary(boundary_2_id) &&
              bulk_el_pt->node_pt(i_node)->is_on_boundary(boundary_1_id))
          {
            // set the output arguments,
            element_pt = dynamic_cast<PERTURBED_ELEMENT*>(bulk_el_pt);
            face_index =
              this->fluid_mesh_pt()->face_index_at_boundary(boundary_1_id, e);

            // Return to exit both loops and end function
            return;
          }
        }
      }
      // If not found, issue warning and return anyway
      oomph_info << "Warning: No corner node found!" << std::endl;
    }

    /// Setup the mesh interactions
    void setup_mesh_interaction()
    {
      SCALING_ELEMENT* singular_el_pt = 0;

      // Loop over the augmented bulk elements
      unsigned n_aug_bulk = Augmented_bulk_element_number.size();
      for (unsigned e = 0; e < n_aug_bulk; e++)
      {
        // Augment elements
        // Upcast from GeneralisedElement to the present element
        PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>(
          this->fluid_mesh_pt()->element_pt(Augmented_bulk_element_number[e]));

        singular_el_pt = dynamic_cast<SCALING_ELEMENT*>(
          Singularity_scaling_mesh_pt->element_pt(0));

        // Set the pointer to the element that determines the amplitude
        // of the singular fct
        el_pt->add_c_equation_element_pt(singular_el_pt);

        singular_el_pt = dynamic_cast<SCALING_ELEMENT*>(
          Singularity_scaling_mesh_pt->element_pt(1));

        // Set the pointer to the element that determines the amplitude
        // of the singular fct
        el_pt->add_c_equation_element_pt(singular_el_pt);
      }
    }
  };
}; // namespace oomph
#endif
