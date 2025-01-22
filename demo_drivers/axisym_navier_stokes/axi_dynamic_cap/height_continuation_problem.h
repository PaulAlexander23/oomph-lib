#ifndef HEIGHT_CONTINUATION_PROBLEM_HEADER
#define HEIGHT_CONTINUATION_PROBLEM_HEADER

#include <fstream>

#include "generic.h"
#include "height_element.h"
#include "temp.h"

namespace oomph
{
  /// A class for the height continuation problem.
  ///
  /// This class is a problem that is used to solve height continuation.
  class HeightContinuationProblem : public Problem, public Temp
  {
  private:
    /// Storage for the height element mesh
    Mesh* Height_mesh_pt;

    /// Storage for the traded height step data
    Data* Traded_data_pt;

  public:
    /// Typedef for the height element
    typedef HeightElement HEIGHT_ELEMENT;

    /// Constructor
    /// Initalises the height element mesh and does nothing else.
    HeightContinuationProblem() : Height_mesh_pt(new Mesh), Traded_data_pt(0)
    {
      add_sub_mesh(Height_mesh_pt);
      build_global_mesh();
    }

    /// Create the height elements using the inner corner and contact line nodes
    /// and trading for the height step data.
    void create_height_elements(SolidNode*& Inner_corner_solid_node_pt,
                                SolidNode*& Contact_line_solid_node_pt,
                                Data*& traded_data_pt)
    {
      // Set the height step data
      Traded_data_pt = traded_data_pt;

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
      // Unset the height step data
      Traded_data_pt = 0;

      // Delete the elements
      unsigned n_element = Height_mesh_pt->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        delete Height_mesh_pt->element_pt(e);
      }
      // Now flush the storage
      Height_mesh_pt->flush_element_and_node_storage();
    }

    /// Takes a continuation step using a fixed height drop and solving for the
    /// parameter.
    double height_step_solve(double ds)
    {
      // Pin the height
      dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
        ->pin_height();
      // Unpin the parameter
      Traded_data_pt->unpin(0);

      // Setup the problem
      oomph_info << "Number of unknowns: " << assign_eqn_numbers() << std::endl;

      // Set the step height
      dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
        ->step_height(ds);

      // Debug the Jacobian
      // this->debug_jacobian();

      // Solve the problem
      this->steady_newton_solve();

      // Restore the problem to its original state
      dynamic_cast<HEIGHT_ELEMENT*>(Height_mesh_pt->element_pt(0))
        ->unpin_height();
      Traded_data_pt->pin(0);

      // Adapt the step size
      ds = adapt_step_size(ds);

      // Return the step size
      return ds;
    }

    /// Adapt the step size based on the number of Newton iterations taken.
    double adapt_step_size(double ds)
    {
      if (this->Nnewton_iter_taken < this->Desired_newton_iterations_ds)
      {
        ds *= 1.5;
      }
      else if (this->Nnewton_iter_taken > this->Desired_newton_iterations_ds)
      {
        ds *= 2.0 / 3.0;
      }
      return ds;
    }

    void set_nnewton_iter_taken(const unsigned& n)
    {
      this->Nnewton_iter_taken = n;
    }

    /// Document the solution
    void doc_solution()
    {
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
