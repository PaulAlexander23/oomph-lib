#ifndef FULL_CONTINUATION_PROBLEM_HEADER
#define FULL_CONTINUATION_PROBLEM_HEADER

#include "singular_axisym_dynamic_cap_problem.h"
#include "height_element.h"
#include "parameters.h"

namespace oomph
{
  template<class ELEMENT, class TIMESTEPPER>
  class FullContinuationProblem
    : public virtual SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>
  {
  private:
    /// Storage for the height element mesh
    Mesh* Height_mesh_pt;

    /// Storage for the traded height step data
    Data* Traded_data_pt;

    /// Flag to indicate if the continuation parameter is set
    bool Is_continuation_parameter_set;

  public:
    /// Typedef for the height element
    typedef HeightElement HEIGHT_ELEMENT;

    /// Constructor
    /// Calls the base constructor and initalises the height element mesh and
    /// traded data.
    FullContinuationProblem(Params* const& params)
      : SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>(params),
        Height_mesh_pt(new Mesh),
        Traded_data_pt(new Data(1)),
        Is_continuation_parameter_set(false)
    {
      this->add_sub_mesh(Height_mesh_pt);
      this->add_global_data(Traded_data_pt);
      this->rebuild_global_mesh();

      Traded_data_pt->pin(0);
    }

    /// actions after adapt
    /// Calls the base actions after adapt and creates the height elements.
    void actions_after_adapt()
    {
      // Call the base actions after adapt
      SingularAxisymDynamicCapProblem<ELEMENT,
                                      TIMESTEPPER>::actions_after_adapt();

      // Create the height elements
      this->create_height_elements(this->inner_corner_solid_node_pt(),
                                   this->contact_line_node_pt());

      // If the continuation parameter is set, pin the height and unpin the
      // parameter
      if (Is_continuation_parameter_set)
      {
        // Pin the height
        dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
          ->pin_height();
        // Unpin the parameter
        Traded_data_pt->unpin(0);


        // Loop over bulk elements and add Reynolds Inverse Froude number as
        // external data
        unsigned n_element = this->bulk_mesh_pt()->nelement();
        for (unsigned e = 0; e < n_element; e++)
        {
          this->bulk_mesh_pt()->element_pt(e)->add_external_data(
            Traded_data_pt);
        }
        // Loop over slip surface elements and add wall velocity as external
        // data
        n_element = this->slip_surface_mesh_pt()->nelement();
        for (unsigned e = 0; e < n_element; e++)
        {
          this->slip_surface_mesh_pt()->element_pt(e)->add_external_data(
            Traded_data_pt);
        }
      }

      // Setup the problem
      this->rebuild_global_mesh();
      oomph_info << "Number of unknowns: " << this->assign_eqn_numbers()
                 << std::endl;
    }

    /// actions before adapt
    /// Calls the base actions before adapt and deletes the height elements.
    void actions_before_adapt()
    {
      // Delete the height elements before the mesh is adapted
      delete_height_elements();

      // Pin the parameter
      Traded_data_pt->pin(0);

      // Call the base actions before adapt, which will also setup the problem
      SingularAxisymDynamicCapProblem<ELEMENT,
                                      TIMESTEPPER>::actions_before_adapt();
    }

    void actions_after_newton_step()
    {
      //this->debug_jacobian();
      //throw OomphLibError(
      //  "Debugging Jacobian", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Set the continuation parameter
    void set_continuation_parameter(double*& parameter_pt)
    {
      // Set the continuation parameter
      Traded_data_pt->set_value(0, *parameter_pt);

      // Delete the original parameter storage
      delete parameter_pt;

      // Use the data storage
      parameter_pt = Traded_data_pt->value_pt(0);

      // set the flag
      Is_continuation_parameter_set = true;
    }

    /// Unset the continuation parameter
    void unset_continuation_parameter(double*& parameter_pt)
    {
      // Create new storage for the parameter and set it to the traded data
      parameter_pt = new double(Traded_data_pt->value(0));

      // Unset the continuation parameter
      Traded_data_pt->set_value(0, 0.0);

      // set the flag
      Is_continuation_parameter_set = false;
    }

    /// Create the height elements using the inner corner and contact line nodes
    /// and trading for the height step data.
    void create_height_elements(SolidNode* Inner_corner_solid_node_pt,
                                SolidNode* Contact_line_solid_node_pt)
    {
      // Create the height element
      HEIGHT_ELEMENT* height_el_pt = new HEIGHT_ELEMENT(
        Inner_corner_solid_node_pt, Contact_line_solid_node_pt);

      // Set the height step data
      height_el_pt->set_parameter_data_pt(Traded_data_pt);

      // Add the height element to the mesh
      Height_mesh_pt->add_element_pt(height_el_pt);
    }

    /// Delete the height elements and unset the height step data.
    void delete_height_elements()
    {
      // Delete the elements
      unsigned n_element = Height_mesh_pt->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        delete Height_mesh_pt->element_pt(e);
      }
      // Now flush the storage
      Height_mesh_pt->flush_element_and_node_storage();
    }

    ///// Takes a continuation step using a fixed height drop and solving for
    /// the
    ///// parameter.
    // double height_step_solve(double ds)
    //{
    //   // Pin the height
    //   dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
    //     ->pin_height();
    //   // Unpin the parameter
    //   Traded_data_pt->unpin(0);

    //  // Setup the problem
    //  oomph_info << "Number of unknowns: " << this->assign_eqn_numbers()
    //             << std::endl;

    //  // Set the step height
    //  dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
    //    ->step_height(ds);

    //  // Debug the Jacobian
    //  // this->debug_jacobian();

    //  // Solve the problem
    //  this->steady_newton_solve();

    //  // Restore the problem to its original state
    //  dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
    //    ->unpin_height();
    //  Traded_data_pt->pin(0);

    //  // Adapt the step size
    //  ds = adapt_step_size(ds);

    //  // Return the step size
    //  return ds;
    //}

    ///// Adapt the step size based on the number of Newton iterations taken.
    // double adapt_step_size(double ds)
    //{
    //   if (this->Nnewton_iter_taken < this->Desired_newton_iterations_ds)
    //   {
    //     ds *= 1.5;
    //   }
    //   else if (this->Nnewton_iter_taken > this->Desired_newton_iterations_ds)
    //   {
    //     ds *= 2.0 / 3.0;
    //   }
    //   return ds;
    // }

    // void set_nnewton_iter_taken(const unsigned& n)
    //{
    //   this->Nnewton_iter_taken = n;
    // }

    void doc_solution()
    {
      SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>::doc_solution();

      // Output the height element
      std::ofstream output_file(this->doc_info().directory() + "/height.dat");
      for (unsigned e = 0; e < Height_mesh_pt->nelement(); e++)
      {
        dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(e))
          ->output(output_file);
      }
      output_file.close();
    }
  };
}; // namespace oomph

#endif
