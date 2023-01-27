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
#include "interface_elements.h"
#include "../generic/integral.h"


namespace oomph
{
  //=========================================================================
  /// Add contribution to element's residual vector and Jacobian
  //=========================================================================
  template<class ELEMENT>
  void FluidInterfaceBoundingElement<ELEMENT>::
    fill_in_generic_residual_contribution_interface_boundary(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Let's get the info from the parent
    ELEMENT* parent_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsids(n_node, dim());
    Vector<double> s_local(dim());

    // Find the dimension of the problem
    unsigned spatial_dim = this->nodal_dimension();

    // Set up the identity matrix
    DenseMatrix<double> identity(spatial_dim, spatial_dim, 0.0);
    for (unsigned i = 0; i < spatial_dim; i++)
    {
      identity(i, i) = 1.0;
    }

    // Storage for the coordinate
    Vector<double> x(spatial_dim);
    Vector<double> dx_dt(spatial_dim);
    Vector<double> u(spatial_dim);

    // Outer unit normal to the wall and surface
    Vector<double> wall_normal(spatial_dim);
    Vector<double> surface_normal(spatial_dim);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Jacobian of mapping
      double J = J_eulerian_at_knot(ipt);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate the shape functions
      shape_at_knot(ipt, psi);

      // Find the dimension of the element
      const unsigned el_dim = dim();
      // Storage for the local coordinates of the integration point
      Vector<double> s_local(el_dim);
      // Set the local coordinate
      for (unsigned i = 0; i < el_dim; i++)
      {
        s_local[i] = integral_pt()->knot(ipt, i);
      }

      // Get the x coordinate
      this->interpolated_x(s_local, x);

      // Get the dx/dt of the coordinate
      const unsigned t_deriv = 1;
      this->interpolated_dxdt(s_local, t_deriv, dx_dt);

      // Assemble interpolated values
      // Loop over the nodes
      double lambda = 0.0;
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          u[i] += nst_u(n, i) * psi(n);
        }
        lambda += lagrange_multiplier(n) * psi(n);
      }
      // Get the unit normal to the wall
      wall_unit_normal(x, wall_normal);

      // Find the local coordinates in the parent
      // Find the dimension of the parent
      unsigned n_dim = parent_pt->dim();
      Vector<double> s_parent(n_dim);
      this->get_local_coordinate_in_bulk(s_local, s_parent);

      // Just get the outer unit normal
      parent_pt->outer_unit_normal(s_parent, surface_normal);

      double contact_angle_local = get_contact_angle(x, dx_dt);

      // Get the value of sigma from the parent
      double sigma_local = parent_pt->sigma(s_parent);

      // Get the capillary number
      double ca_local = ca();

      // Get the wall tangent vector
      Vector<double> wall_tangent(spatial_dim);
      // Rotate the normal vector clockwise 90 degrees
      wall_tangent[0] = wall_normal[1];
      wall_tangent[1] = -wall_normal[0];

      Vector<double> project_to_interface_normal_on_wall(spatial_dim, 0.0);
      for (unsigned d = 0; d < spatial_dim; d++)
      {
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          project_to_interface_normal_on_wall[d] +=
            (identity(d, i) - wall_normal[d] * wall_normal[i]) *
            surface_normal[d];
        }
      }

      // Loop over the nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // Add the tension component to the momentum equations
        // Used in both the strong and weak forms
        for (unsigned i = 0; i < spatial_dim; i++)
        {
          int local_eqn = nst_momentum_local_eqn(n, i);
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -= (sigma_local / ca_local) *
                                    (sin(contact_angle_local) * wall_normal[i] +
                                     cos(contact_angle_local) *
                                       project_to_interface_normal_on_wall[i]);
          }
        }

        // STRONG contact angle imposition
        if (Contact_angle_flag[n] == 1)
        {
          /// Strong form of the contact angle constraint by hijacking of the
          /// free surface kinematic equation

          // Read out the kinematic equation number
          int local_eqn = this->fsi_kinematic_local_eqn(n);

          // Note that because we have outer unit normals for the free surface
          // and the wall, the cosine of the contact angle is equal to
          // MINUS the dot product computed above
          if (local_eqn >= 0)
          {
            // Find the dot product of the two vectors
            double dot = 0.0;
            for (unsigned i = 0; i < spatial_dim; i++)
            {
              dot += surface_normal[i] * wall_normal[i];
            }

            // Reset the residuals, as we are hijacking.
            // TODO Check that we need to do this.
            if (n == 0)
            {
              residuals[local_eqn] = 0.0;
            }

            residuals[local_eqn] +=
              (cos(contact_angle_local) + dot) * psi(n) * W;
          }
        }
        // WEAK contact angle imposition
        else if (Contact_angle_flag[n] == 2)
        {
          int local_eqn = wall_bounded_kinematic_local_eqn(n);
          if (local_eqn >= 0)
          {
            // Point kinematic equation
            double st_local = parent_pt->st();

            // Wall-bounded kinematic equation
            for (unsigned k = 0; k < spatial_dim; k++)
            {
              residuals[local_eqn] += (u[k] - st_local * dx_dt[k]) *
                                      project_to_interface_normal_on_wall[k] *
                                      psi(n) * W;

              // Lagrange multiplier contribution to momentum equations
              local_eqn = this->nst_momentum_local_eqn(n, k);
              if (local_eqn >= 0)
              {
                residuals[local_eqn] +=
                  lambda * project_to_interface_normal_on_wall[k] * psi(n) * W;
              }

              // Lagrange multiplier contribution to solid displacement
              // equations
              local_eqn = this->fsi_kinematic_local_eqn(n);
              if (local_eqn >= 0)
              {
                double residual_contribution =
                  -lambda * st_local * project_to_interface_normal_on_wall[k] *
                  psi(n) * W;
                // TODO Check if we need further contributions from
                // project_to_interface_normal_on_wall
                // If we are timestepping,...
                const unsigned time_derivative = 1;
                const unsigned weight_number = 0;
                const double time_stepper_weight =
                  this->node_pt(n)->time_stepper_pt()->weight(time_derivative,
                                                              weight_number);
                if (time_stepper_weight > 0)
                {
                  // ... scale the contribution by the weight.
                  residual_contribution *= time_stepper_weight;
                }

                residuals[local_eqn] += residual_contribution;
              }
            }
          }
        }
        // NOTE: The jacobian entries will be computed automatically
        // by finite differences.

      } // Node loop

      // Now add the additional contributions
      add_additional_residual_contributions_interface_boundary(
        residuals, jacobian, flag, psi, dpsids, surface_normal, W);
    } // Integration points loop
  }


  /// //////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////

  //====================================================================
  /// Specialise the surface derivatives for the line interface case
  //===================================================================
  double LineDerivatives::compute_surface_derivatives(
    const Shape& psi,
    const DShape& dpsids,
    const DenseMatrix<double>& interpolated_t,
    const Vector<double>& interpolated_x,
    DShape& dpsidS,
    DShape& dpsidS_div)
  {
    const unsigned n_shape = psi.nindex1();
    const unsigned n_dim = 2;

    // Calculate the only entry of the surface
    // metric tensor (square length of the tangent vector)
    double a11 = interpolated_t(0, 0) * interpolated_t(0, 0) +
                 interpolated_t(0, 1) * interpolated_t(0, 1);

    // Now set the derivative terms of the shape functions
    for (unsigned l = 0; l < n_shape; l++)
    {
      for (unsigned i = 0; i < n_dim; i++)
      {
        dpsidS(l, i) = dpsids(l, 0) * interpolated_t(0, i) / a11;
      }
    }

    // The surface divergence is the same as the surface
    // gradient operator
    dpsidS_div = dpsidS;

    // Return the jacobian
    return sqrt(a11);
  }


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////

  //====================================================================
  /// Specialise the surface derivatives for the axisymmetric interface case
  //===================================================================
  double AxisymmetricDerivatives::compute_surface_derivatives(
    const Shape& psi,
    const DShape& dpsids,
    const DenseMatrix<double>& interpolated_t,
    const Vector<double>& interpolated_x,
    DShape& dpsidS,
    DShape& dpsidS_div)
  {
    // Initially the same as the 2D case
    const unsigned n_shape = psi.nindex1();
    const unsigned n_dim = 2;

    // Calculate the only entry of the surface
    // metric tensor (square length of the tangent vector)
    double a11 = interpolated_t(0, 0) * interpolated_t(0, 0) +
                 interpolated_t(0, 1) * interpolated_t(0, 1);

    // Now set the derivative terms of the shape functions
    for (unsigned l = 0; l < n_shape; l++)
    {
      for (unsigned i = 0; i < n_dim; i++)
      {
        dpsidS(l, i) = dpsids(l, 0) * interpolated_t(0, i) / a11;
        // Set the standard components of the divergence
        dpsidS_div(l, i) = dpsidS(l, i);
      }
    }

    const double r = interpolated_x[0];

    // The surface divergence is different because we need
    // to take account of variation of the base vectors
    for (unsigned l = 0; l < n_shape; l++)
    {
      dpsidS_div(l, 0) += psi(l) / r;
    }

    // Return the jacobian, needs to be multiplied by r
    return r * sqrt(a11);
  }


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////

  //====================================================================
  /// Specialise the surface derivatives for 2D surface case
  //===================================================================
  double SurfaceDerivatives::compute_surface_derivatives(
    const Shape& psi,
    const DShape& dpsids,
    const DenseMatrix<double>& interpolated_t,
    const Vector<double>& interpolated_x,
    DShape& dpsidS,
    DShape& dpsidS_div)
  {
    const unsigned n_shape = psi.nindex1();
    const unsigned n_dim = 3;

    // Calculate the local metric tensor
    // The dot product of the two tangent vectors
    double amet[2][2];
    for (unsigned al = 0; al < 2; al++)
    {
      for (unsigned be = 0; be < 2; be++)
      {
        // Initialise to zero
        amet[al][be] = 0.0;
        // Add the dot product contributions
        for (unsigned i = 0; i < n_dim; i++)
        {
          amet[al][be] += interpolated_t(al, i) * interpolated_t(be, i);
        }
      }
    }

    // Work out the determinant
    double det_a = amet[0][0] * amet[1][1] - amet[0][1] * amet[1][0];
    // Also work out the contravariant metric
    double aup[2][2];
    aup[0][0] = amet[1][1] / det_a;
    aup[0][1] = -amet[0][1] / det_a;
    aup[1][0] = -amet[1][0] / det_a;
    aup[1][1] = amet[0][0] / det_a;


    // Now construct the surface gradient terms
    for (unsigned l = 0; l < n_shape; l++)
    {
      // Do some pre-computation
      const double dpsi_temp[2] = {
        aup[0][0] * dpsids(l, 0) + aup[0][1] * dpsids(l, 1),
        aup[1][0] * dpsids(l, 0) + aup[1][1] * dpsids(l, 1)};

      for (unsigned i = 0; i < n_dim; i++)
      {
        dpsidS(l, i) = dpsi_temp[0] * interpolated_t(0, i) +
                       dpsi_temp[1] * interpolated_t(1, i);
      }
    }

    // The divergence operator is the same
    dpsidS_div = dpsidS;

    // Return the jacobian
    return sqrt(det_a);
  }

} // namespace oomph
