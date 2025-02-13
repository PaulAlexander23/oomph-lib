#ifndef SINGULAR_OVERLAYING_MY_LINEAR_ELEMENT_HEADER
#define SINGULAR_OVERLAYING_MY_LINEAR_ELEMENT_HEADER

#include "overlaying_my_linear_element.h"

namespace oomph
{
  template<class BASE_ELEMENT>
  class SingularOverlayingMyLinearElement
    : public virtual OverlayingMyLinearElement<BASE_ELEMENT>
  {
  private:
    /// Flag to indicate if the element has been augmented
    bool IsAugmented;

    /// Flag to indicate if the Jacobian is computed using finite differences
    bool IsJacobianFD;

    /// Vector of pointers to SingularNavierStokesSolutionElement objects
    Vector<SingularNavierStokesSolutionElement<
      OverlayingMyLinearElement<BASE_ELEMENT>>*>
      C_equation_elements_pt;

  public:
    SingularOverlayingMyLinearElement()
      : OverlayingMyLinearElement<BASE_ELEMENT>(),
        IsAugmented(false),
        IsJacobianFD(false)
    {
    }

    bool is_augmented() const
    {
      return IsAugmented;
    }

    // Check if the element is using finite differences for the Jacobian
    bool is_using_fd_jacobian()
    {
      return IsJacobianFD;
    }

    void augment()
    {
      // Loop over the nodes and add extra data
      const unsigned n_node = this->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        const unsigned original_n_value = this->node_pt(n)->nvalue();
        // We need an additional n_u_lin_axi_nst values for the total velocity
        // equations
        unsigned desired_n_value =
          this->n_r_lin_el() + 2 * this->n_u_lin_axi_nst();

        // If this is a pressure node then we need to add additional pressure
        // values.
        if (this->is_pressure_node(n))
        {
          desired_n_value += 2 * this->n_p_lin_axi_nst();
        }

        // If we need to
        if (original_n_value != desired_n_value)
        {
          // resize
          this->node_pt(n)->resize(desired_n_value);
        }

        // Temporarily pin the additional values
        for (unsigned i = original_n_value; i < desired_n_value; i++)
        {
          this->node_pt(n)->pin(i);
        }
      }

      IsAugmented = true;
    }

    virtual inline unsigned u_index_lin_axi_nst_fe(const unsigned& n,
                                                   const unsigned& i)
    {
      unsigned index = this->n_r_lin_el() + this->n_u_lin_axi_nst() + i;
      if (this->is_pressure_node(n))
      {
        index += this->n_p_lin_axi_nst();
      }
      return index;
    }

    virtual inline unsigned p_index_lin_axi_nst_fe(const unsigned& n,
                                                   const unsigned& i)
    {
      return this->n_r_lin_el() + 2 * this->n_u_lin_axi_nst() + this->n_p_lin_axi_nst() + i;
    }

    void add_c_equation_element_pt(
      SingularNavierStokesSolutionElement<
        OverlayingMyLinearElement<BASE_ELEMENT>>* c_pt)
    {
      // Add the element
      C_equation_elements_pt.push_back(c_pt);

      // Add the additional unknown of this object as external data in the
      // Navier-Stokes element
      const bool use_fd = false;
      this->add_external_data(c_pt->internal_data_pt(0), use_fd);
    }

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      // Get the contribution from the underlying wrapped element first
      OverlayingMyLinearElement<
        BASE_ELEMENT>::fill_in_contribution_to_residuals(residuals);

      if (this->is_augmented())
      {
        this->fill_in_generic_residual_contribution_wrapped_axi_nst(
          residuals, GeneralisedElement::Dummy_matrix, 0);
      }
    }

    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Use finite differences
      if (this->is_using_fd_jacobian())
      {
        FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      }
      // Otherwise use analytic contributions
      else
      {
        // Call the base fill_in_contribution_to_jacobian function
        OverlayingMyLinearElement<
          BASE_ELEMENT>::fill_in_contribution_to_jacobian(residuals, jacobian);
        // Then call the singular Navier-Stokes element's
        // fill_in_contribution_to_jacobian function
        if (this->is_augmented())
        {
          this->fill_in_generic_residual_contribution_wrapped_axi_nst(
            residuals, jacobian, 1);
        }
      }
    }

    /// Fill in the additional contributions to the momentum equations and
    /// implement the total velocity equations
    void fill_in_generic_residual_contribution_wrapped_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
    }

    /// Return the i-th component of the FE interpolated velocity
    /// u[i] at local coordinate s
    double interpolated_u_lin_axi_nst_fe(const Vector<double>& s,
                                         const unsigned& i)
    {
      // Determine number of nodes in the element
      const unsigned n_node = this->nnode();

      // Provide storage for local shape functions
      Shape psi(n_node);

      // Find values of shape functions
      this->shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(l, u_index_lin_axi_nst_fe(l, i)) * psi[l];
      }

      return (interpolated_u);
    }

    /// Return the i-th component of the FE interpolated pressure
    /// p[i] at local coordinate s
    double interpolated_p_lin_axi_nst_fe(const Vector<double>& s,
                                         const unsigned& i)
    {
      // Determine number of pressure nodes in the element
      const unsigned n_pressure_nodes = this->npres_lin_axi_nst();

      // Provide storage for local shape functions
      Shape psi(n_pressure_nodes);

      // Find values of shape functions
      this->pshape_lin_axi_nst(s, psi);

      // Initialise value of p
      double interpolated_p = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_pressure_nodes; l++)
      {
        // N.B. The pure virtual function p_lin_axi_nst(...)
        // automatically calculates the index at which the pressure value
        // is stored, so we don't need to worry about this here
        interpolated_p += this->nodal_value(l, p_index_lin_axi_nst_fe(l, i)) * psi[l];
      }

      return (interpolated_p);
    }

    double interpolated_u_lin_axi_nst_bar(const Vector<double>& s,
                                          const unsigned& i)
    {
      return 0.0;
    }

    double interpolated_p_lin_axi_nst_bar(const Vector<double>& s,
                                          const unsigned& i)
    {
      return 0.0;
    }


    /// Output function in tecplot format:
    /// r, z,
    /// Displacements: R^C, R^S, Z^C, Z^S,
    /// total values: U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
    /// fe correction values: U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S,
    /// singular solution: U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S,
    /// Error, Size
    /// Specified number of plot points in each coordinate direction.
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Provide storage for vector of local coordinates
      Vector<double> s(2);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Determine number of plot points
      const unsigned n_plot_points = this->nplot_points(nplot);

      // Loop over plot points
      for (unsigned iplot = 0; iplot < n_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        // Output global coordinates to file
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }

        // Output perturbations to nodal positions to file
        for (unsigned i = 0; i < 4; i++)
        {
          outfile << this->interpolated_nodal_position_perturbation_lin_axi_nst(s, i)
                  << " ";
        }

        //  Output velocities to file
        for (unsigned i = 0; i < 6; i++)
        {
          outfile << this->interpolated_u_lin_axi_nst(s, i) << " ";
        }

        // Output pressure to file
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << this->interpolated_p_lin_axi_nst(s, i) << " ";
        }

        //  Output velocities to file
        for (unsigned i = 0; i < 6; i++)
        {
          outfile << interpolated_u_lin_axi_nst_fe(s, i) << " ";
        }

        // Output pressure to file
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << interpolated_p_lin_axi_nst_fe(s, i) << " ";
        }

        //  Output velocities to file
        for (unsigned i = 0; i < 6; i++)
        {
          outfile << interpolated_u_lin_axi_nst_bar(s, i) << " ";
        }

        // Output pressure to file
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << interpolated_p_lin_axi_nst_bar(s, i) << " ";
        }

        // Error
        outfile << this->error() << " ";

        // Size
        outfile << this->size() << " ";

        outfile << std::endl;
      }
      outfile << std::endl;

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);

    } // End of output
  };

  template<class BASE_ELEMENT>
  class FaceGeometry<SingularOverlayingMyLinearElement<BASE_ELEMENT>>
    : public TElement<1, 3>
  {
  };
}; // namespace oomph

#endif
