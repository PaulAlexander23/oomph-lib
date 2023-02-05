// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Non-inline functions for fluid free surface elements

// OOMPH-LIB headers
#include "interface_contact_elements.h"


namespace oomph
{
  //=========================================================================
  /// Set a pointer to the desired contact angle. Optional boolean
  /// (defaults to true)
  /// chooses strong imposition via hijacking (true) or weak imposition
  /// via addition to momentum equation (false). The default strong imposition
  /// is appropriate for static contact problems.
  //=========================================================================
  void FluidInterfaceBoundingElement::set_contact_angle(double* const& angle_pt,
                                                        const bool& strong)
  {
    // Set the pointer to the contact angle
    Contact_angle_pt = angle_pt;

    // Set the contact angle fct flag to use the fixed value
    Contact_angle_fct_flag = 0;

    // If we are hijacking the kinematic condition (the default)
    // to do the strong (pointwise form of the contact-angle condition)
    if (strong)
    {
      // Remember what we're doing
      Contact_angle_flag = 1;

      // Hijack the bulk element residuals
      dynamic_cast<FluidInterfaceElement*>(bulk_element_pt())
        ->hijack_kinematic_conditions(Bulk_node_number);
    }
    // Otherwise, we'll impose it weakly via the momentum equations.
    // This will require that the appropriate velocity node is unpinned,
    // which is why this is a bad choice for static contact problems in which
    // there is a no-slip condition on the wall. In that case, the momentum
    // equation is never assembled and so the contact angle condition is not
    // applied unless we use the strong version above.
    else
    {
      Contact_angle_flag = 2;
    }
  }

  //=========================================================================
  /// Set a pointer to the desired contact angle function. Optional boolean
  /// (defaults to true)
  /// chooses strong imposition via hijacking (true) or weak imposition
  /// via addition to momentum equation (false). The default strong imposition
  /// is appropriate for static contact problems.
  //=========================================================================
  void FluidInterfaceBoundingElement::set_contact_angle_fct(
    ContactAngleFctPt const& angle_fct_pt, const bool& strong)
  {
    // Set the pointer to the contact angle
    Contact_angle_fct_pt = angle_fct_pt;

    // Set the contact angle fct flag to use the fixed value
    Contact_angle_fct_flag = 1;

    // If we are hijacking the kinematic condition (the default)
    // to do the strong (pointwise form of the contact-angle condition)
    if (strong)
    {
      // Remember what we're doing
      Contact_angle_flag = 1;

      // Hijack the bulk element residuals
      dynamic_cast<FluidInterfaceElement*>(bulk_element_pt())
        ->hijack_kinematic_conditions(Bulk_node_number);
    }
    // Otherwise, we'll impose it weakly via the momentum equations.
    // This will require that the appropriate velocity node is unpinned,
    // which is why this is a bad choice for static contact problems in which
    // there is a no-slip condition on the wall. In that case, the momentum
    // equation is never assembled and so the contact angle condition is not
    // applied unless we use the strong version above.
    else
    {
      Contact_angle_flag = 2;
    }
  }

  //=========================================================================
  /// Add contribution to element's residual vector and Jacobian
  //=========================================================================
  void PointFluidInterfaceBoundingElement::
    fill_in_generic_residual_contribution_interface_boundary(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Let's get the info from the parent
    FiniteElement* parent_pt = bulk_element_pt();

    // Find the dimension of the problem
    unsigned spatial_dim = this->nodal_dimension();

    // Outer unit normal to the wall
    Vector<double> wall_normal(spatial_dim);

    // Outer unit normal to the free surface
    Vector<double> unit_normal(spatial_dim);

    // Storage for the coordinate
    Vector<double> x(spatial_dim);
    double v = 0.0;

    // Storage for the coordinate time derivative
    Vector<double> dx_dt(spatial_dim);

    // Find the dimension of the parent
    unsigned n_dim = parent_pt->dim();

    // Dummy local coordinate, of size zero
    Vector<double> s_local(0);

    // Get the x coordinate
    this->interpolated_x(s_local, x);

    // Get the dx/dt of the coordinate
    const unsigned t_deriv = 1;
    this->interpolated_dxdt(s_local, t_deriv, dx_dt);

    // Get the unit normal to the wall
    wall_unit_normal(x, wall_normal);

    // Find the local coordinates in the parent
    Vector<double> s_parent(n_dim);
    this->get_local_coordinate_in_bulk(s_local, s_parent);

    // Just get the outer unit normal
    dynamic_cast<FaceElement*>(parent_pt)->outer_unit_normal(s_parent,
                                                             unit_normal);

    v = dynamic_cast<FluidInterfaceElement*>(this->bulk_element_pt())
          ->interpolated_u(s_parent, 1);


    // Find the dot product of the two vectors
    double dot = 0.0;
    for (unsigned i = 0; i < spatial_dim; i++)
    {
      dot += unit_normal[i] * wall_normal[i];
    }

    // Get the value of sigma from the parent
    double sigma_local =
      dynamic_cast<FluidInterfaceElement*>(parent_pt)->sigma(s_parent);

    // Get the capillary number
    double ca_local = ca();

    double contact_angle_local = 0.0;
    if (Contact_angle_fct_flag)
    {
      contact_angle(dx_dt, contact_angle_local);
    }
    else
    {
      contact_angle(contact_angle_local);
    }


    // Add the tension component to the momentum equations
    // Used in both the strong and weak forms

    // Get the wall tangent vector
    Vector<double> wall_tangent(spatial_dim);
    // Rotate the normal vector clockwise 90 degrees
    wall_tangent[0] = wall_normal[1];
    wall_tangent[1] = -wall_normal[0];

    //// Get the wall tangent vector
    // Vector<double> unit_tangent(spatial_dim);
    //// Rotate the normal vector clockwise 90 degrees
    // unit_tangent[0] = unit_normal[1];
    // unit_tangent[1] = -unit_normal[0];

    //// Add the dynamic boundary condition contribution to the momentum
    /// equations
    // for (unsigned i = 0; i < 2; i++)
    //{
    //  int local_eqn = nodal_local_eqn(0,
    //  this->U_index_interface_boundary[0][i]); if (local_eqn >= 0)
    //  {
    //    residuals[local_eqn] += (sigma_local / ca_local) * unit_tangent[i];
    //  }
    //}

    // Just add the appropriate contribution to the momentum equations
    for (unsigned i = 0; i < 2; i++)
    {
      int local_eqn =
        nodal_local_eqn(0, this->U_index_interface_boundary[0][i]);
      if (local_eqn >= 0)
      {
        residuals[local_eqn] -= (sigma_local / ca_local) *
                                (sin(contact_angle_local) * wall_normal[i] +
                                 cos(contact_angle_local) * wall_tangent[i]);
      }
    }


    switch (Contact_angle_flag)
    {
        /// Strong form of the contact angle constraint by hijacking of the
        /// kinematic equation
      case 1:
      {
        // If we are imposing the contact angle strongly (by hijacking)
        // overwrite the kinematic equation

        // Read out the kinematic equation number
        int local_eqn = kinematic_local_eqn(0);

        // Note that because we have outer unit normals for the free surface
        // and the wall, the cosine of the contact angle is equal to
        // MINUS the dot product computed above
        if (local_eqn >= 0)
        {
          residuals[local_eqn] = cos(contact_angle_local) + dot;
        }

        break;
      }
        /// Weak form
      case 2:
      {
        // The contribution to the momentum equations suffices, so we don't
        // need to do anything else for the weak form

        break;
      }
        /// Default - Shouldn't get here
      default:
        break;
    }
    // NOTE: The jacobian entries will be computed automatically
    // by finite differences.


    // Point kinematic equation
    if (Kinematic_lagrange_index >= 0)
    {
      double lambda = internal_data_pt(Kinematic_lagrange_index)->value(0);
      int local_eqn = internal_local_eqn(Kinematic_lagrange_index, 0);
      if (local_eqn >= 0)
      {
        double st_local = dynamic_cast<FluidInterfaceElement*>(parent_pt)->st();
        // Constraint
        residuals[local_eqn] = v - st_local * dx_dt[1];

        // Lagrange multiplier contribution to momentum equations
        local_eqn =
          this->nodal_local_eqn(0, this->U_index_interface_boundary[0][1]);
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += lambda;
        }

        // Lagrange multiplier contribution to solid displacement equations
        local_eqn = this->kinematic_local_eqn(0);
        if (local_eqn >= 0)
        {
          if (this->node_pt(0)->time_stepper_pt()->weight(1, 0) > 0)
          {
            residuals[local_eqn] -=
              lambda * st_local *
              this->node_pt(0)->time_stepper_pt()->weight(1, 0);
          }
          else
          {
            residuals[local_eqn] -= lambda;
          }
        }
      }
    }


    // Dummy arguments
    Shape psif(1);
    DShape dpsifds(1, 1);
    Vector<double> interpolated_n(1);
    double W = 0.0;

    // Now add the additional contributions
    add_additional_residual_contributions_interface_boundary(
      residuals, jacobian, flag, psif, dpsifds, interpolated_n, W);
  }


  //=========================================================================
  /// Add contribution to element's residual vector and Jacobian
  //=========================================================================
  void LineFluidInterfaceBoundingElement::
    fill_in_generic_residual_contribution_interface_boundary(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Let's get the info from the parent
    FiniteElement* parent_pt = bulk_element_pt();

    // Find the dimension of the problem
    unsigned spatial_dim = this->nodal_dimension();

    // Outer unit normal to the wall
    Vector<double> wall_normal(spatial_dim);

    // Outer unit normal to the free surface
    Vector<double> unit_normal(spatial_dim);

    // Find the dimension of the parent
    unsigned n_dim = parent_pt->dim();

    // Find the local coordinates in the parent
    Vector<double> s_parent(n_dim);

    // Storage for the shape functions
    unsigned n_node = this->nnode();
    Shape psi(n_node);
    DShape dpsids(n_node, 1);
    Vector<double> s_local(1);

    // Loop over intergration points
    unsigned n_intpt = this->integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ++ipt)
    {
      // Get the local coordinate of the integration point
      s_local[0] = this->integral_pt()->knot(ipt, 0);
      get_local_coordinate_in_bulk(s_local, s_parent);

      // Get the local shape functions
      this->dshape_local(s_local, psi, dpsids);

      // Zero the position
      Vector<double> x(spatial_dim, 0.0);

      // Now construct the position and the tangent
      Vector<double> interpolated_t1(spatial_dim, 0.0);
      for (unsigned n = 0; n < n_node; n++)
      {
        const double psi_local = psi(n);
        const double dpsi_local = dpsids(n, 0);
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          double pos = this->nodal_position(n, i);
          interpolated_t1[i] += pos * dpsi_local;
          x[i] += pos * psi_local;
        }
      }

      // Now we can calculate the Jacobian term
      double t_length = 0.0;
      for (unsigned i = 0; i < spatial_dim; ++i)
      {
        t_length += interpolated_t1[i] * interpolated_t1[i];
      }
      double W = std::sqrt(t_length) * this->integral_pt()->weight(ipt);

      // Storage for the coordinate time derivative
      Vector<double> dx_dt(spatial_dim);

      // Get the dx/dt of the coordinate
      const unsigned t_deriv = 1;
      this->interpolated_dxdt(s_local, t_deriv, dx_dt);

      // Get the contact angle
      double contact_angle_local = 0.0;
      if (Contact_angle_fct_flag)
      {
        contact_angle(dx_dt, contact_angle_local);
      }
      else
      {
        contact_angle(contact_angle_local);
      }

      // Imposition of contact angle in weak form
      if (Contact_angle_flag == 2)
      {
        // Get the outer unit normal of the entire interface
        dynamic_cast<FaceElement*>(parent_pt)->outer_unit_normal(s_parent,
                                                                 unit_normal);

        // Calculate the wall normal
        wall_unit_normal(x, wall_normal);

        // Find the dot product of the two
        double dot = 0.0;
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          dot += unit_normal[i] * wall_normal[i];
        }

        // Find the projection of the outer normal of the surface into the plane
        Vector<double> binorm(spatial_dim);
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          binorm[i] = unit_normal[i] - dot * wall_normal[i];
        }

        // Get the value of sigma from the parent
        const double sigma_local =
          dynamic_cast<FluidInterfaceElement*>(parent_pt)->sigma(s_parent);

        // Get the capillary number
        const double ca_local = ca();


        // Add the contributions to the momentum equation

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < 3; i++)
          {
            // Get the equation number for the momentum equation
            int local_eqn =
              this->nodal_local_eqn(l, this->U_index_interface_boundary[l][i]);

            // If it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Add the surface-tension contribution to the momentum equation
              residuals[local_eqn] +=
                (sigma_local / ca_local) *
                (sin(contact_angle_local) * wall_normal[i] +
                 cos(contact_angle_local) * binorm[i]) *
                psi(l) * W;
            }
          }
        }
      }
      // Otherwise [strong imposition (by hijacking) of contact angle or
      // "no constraint at all"], add the appropriate contribution to
      // the momentum equation
      else
      {
        // Storage for the outer vector
        Vector<double> m(3);

        // Get the outer unit normal of the line
        this->outer_unit_normal(s_local, m);

        // Get the value of sigma from the parent
        const double sigma_local =
          dynamic_cast<FluidInterfaceElement*>(parent_pt)->sigma(s_parent);

        // Get the capillary number
        const double ca_local = ca();

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < 3; i++)
          {
            // Get the equation number for the momentum equation
            int local_eqn =
              this->nodal_local_eqn(l, this->U_index_interface_boundary[l][i]);

            // If it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Add the surface-tension contribution to the momentum equation
              residuals[local_eqn] +=
                m[i] * (sigma_local / ca_local) * psi(l) * W;
            }
          }
        }
      } // End of the line integral terms


      // If we are imposing the contact angle strongly (by hijacking)
      // overwrite the kinematic equation
      if (Contact_angle_flag == 1)
      {
        // Get the outer unit normal of the whole interface
        dynamic_cast<FaceElement*>(parent_pt)->outer_unit_normal(s_parent,
                                                                 unit_normal);

        // Calculate the wall normal
        wall_unit_normal(x, wall_normal);

        // Find the dot product
        double dot = 0.0;
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          dot += unit_normal[i] * wall_normal[i];
        }

        // Loop over the test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Read out the kinematic equation number
          int local_eqn = kinematic_local_eqn(l);

          // Note that because we have outer unit normals for the free surface
          // and the wall, the cosine of the contact angle is equal to
          // MINUS the dot product
          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              (cos(contact_angle_local) + dot) * psi(l) * W;
          }
          // NOTE: The jacobian entries will be computed automatically
          // by finite differences.
        }
      } // End of strong form of contact angle condition

      // Add any additional residual contributions
      add_additional_residual_contributions_interface_boundary(
        residuals, jacobian, flag, psi, dpsids, unit_normal, W);
    }
  }


} // namespace oomph
