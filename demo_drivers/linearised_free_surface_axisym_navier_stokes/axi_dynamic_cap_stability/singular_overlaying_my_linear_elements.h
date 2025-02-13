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
        IsJacobianFD(true)
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
      return this->n_r_lin_el() + 2 * this->n_u_lin_axi_nst() +
             this->n_p_lin_axi_nst() + i;
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

    /// Evaluate sum of all velocity singular fcts
    /// (incl. the amplitude) at Eulerian position x
    double u_bar(const Vector<double>& x, const unsigned& i) const
    {
      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the dimension of the problem
      double sum = 0.0;
      for (unsigned s = 0; s < n_sing; s++)
      {
        Vector<double> u_bar_local = C_equation_elements_pt[s]->u_bar(x);
        sum += u_bar_local[i];
      }
      return sum;
    }

    /// Evaluate gradient of sum of all velocity singular fcts
    /// (incl. the amplitudes) at Eulerian position x: grad[i][j] = du_i/dx_j
    Vector<Vector<double>> grad_u_bar(const Vector<double>& x) const
    {
      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      // Find the dimension of the problem
      unsigned cached_dim = this->dim();
      Vector<Vector<double>> sum(cached_dim);
      for (unsigned i = 0; i < cached_dim; i++)
      {
        sum[i].resize(cached_dim, 0.0);
      }
      for (unsigned s = 0; s < n_sing; s++)
      {
        Vector<Vector<double>> grad_u_bar_local =
          C_equation_elements_pt[s]->grad_u_bar(x);
        for (unsigned i = 0; i < cached_dim; i++)
        {
          for (unsigned j = 0; j < cached_dim; j++)
          {
            sum[i][j] += grad_u_bar_local[i][j];
          }
        }
      }
      return sum;
    }

    /// Evaluate sum of all pressure singular fcts
    /// (incl. the amplitudes) at Eulerian position x
    double p_bar(const Vector<double>& x) const
    {
      // Find the number of singularities
      unsigned n_sing = C_equation_elements_pt.size();

      double sum = 0.0;
      for (unsigned i = 0; i < n_sing; i++)
      {
        sum += C_equation_elements_pt[i]->p_bar(x);
      }
      return sum;
    }

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      // Get the contribution from the underlying wrapped element first

      if (this->is_augmented())
      {
        this->fill_in_generic_residual_contribution_wrapped_axi_nst(
          residuals,
          GeneralisedElement::Dummy_matrix,
          GeneralisedElement::Dummy_matrix,
          0);
        this->fill_in_generic_contribution_to_residuals_linear_elasticity(
          residuals, GeneralisedElement::Dummy_matrix, 0);
      }
      else
      {
        OverlayingMyLinearElement<
          BASE_ELEMENT>::fill_in_contribution_to_residuals(residuals);
      }
    }

    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      this->fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Fill in the additional contributions to the momentum equations and
    /// implement the total velocity equations
    void fill_in_generic_residual_contribution_wrapped_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      const unsigned& flag)
    {
      // Get the time from the first node in the element
      const double time = this->node_pt(0)->time_stepper_pt()->time();

      // Determine number of nodes in the element
      const unsigned n_node = this->nnode();

      // Determine how many pressure values there are associated with
      // a single pressure component
      const unsigned n_pres = this->npres_lin_axi_nst();

      // Get the nodal indices at which the velocity is stored
      unsigned u_nodal_index[6];
      for (unsigned i = 0; i < 6; ++i)
      {
        u_nodal_index[i] = this->u_index_lin_axi_nst(i);
      }

      // Set up memory for the fluid shape and test functions
      // Note that there are two spatial dimensions, r and z, in this problem
      Shape psif(n_node), testf(n_node);
      DShape dpsifds(n_node, 2), dpsifdx(n_node, 2);
      DShape dtestfds(n_node, 2), dtestfdx(n_node, 2);

      // Set up memory for the pressure (p) shape and test functions
      Shape psip(n_pres), testp(n_pres);

      // Determine number of integration points
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Set up memory for the vector to hold local coordinates (two dimensions)
      Vector<double> s(2);

      // Get physical variables from the element
      // (Reynolds number must be multiplied by the density ratio)
      const double scaled_re = this->re() * this->density_ratio();
      const double scaled_re_st = this->re_st() * this->density_ratio();
      const double scaled_re_inv_fr = this->re_invfr() * this->density_ratio();
      const double visc_ratio = this->viscosity_ratio();
      Vector<double> G = this->g();
      const int k = this->azimuthal_mode_number();

      // Integers used to store the local equation and unknown numbers
      int local_eqn = 0, local_unknown = 0;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of the local coordinates s
        for (unsigned i = 0; i < 2; i++)
        {
          s[i] = this->integral_pt()->knot(ipt, i);
        }

        // Get the integral weight
        const double w = this->integral_pt()->weight(ipt);

        // Calculate the derivatives of the fluid shape functions w.r.t. the
        // local coordinates
        this->dshape_local_at_knot(ipt, psif, dpsifds);

        // Set the derivatives of the test functions w.r.t. local coords
        // equal to that of the shape functions
        dtestfds = dpsifds;

        // Calculate the fluid shape and test functions, and their derivatives
        // w.r.t. the global coordinates
        const double Jbar = this->dshape_and_dtest_eulerian_at_knot_lin_axi_nst(
          ipt, psif, dpsifdx, testf, dtestfdx);

        // Calculate the pressure shape and test functions
        this->pshape_lin_axi_nst(s, psip, testp);

        // Allocate storage for the position and the derivative of the
        // mesh positions w.r.t. time
        Vector<double> interpolated_x(2, 0.0);
        Vector<double> mesh_velocity(2, 0.0);

        // Allocate storage for the derivatives of the unperturbed positions
        // w.r.t. local coordinates (s_1 and s_2)
        DenseMatrix<double> interpolated_dxbar_ds(2, 2, 0.0);

        // Allocate storage for the perturbed position xhat
        Vector<double> interpolated_xhat(4, 0.0);

        // Allocate storage for the derivative of the perturbed position
        // w.r.t. time
        Vector<double> dxhat_dt(4, 0.0);

        // Allocate storage for the derivatives of the perturbed positions
        // w.r.t. local coordinates (s_1 and s_2)
        DenseMatrix<double> interpolated_dxhat_ds(4, 2, 0.0);

        // Allocate storage for the velocity components (six of these)
        // and their derivatives w.r.t. time
        Vector<double> interpolated_u(6, 0.0);
        Vector<double> dudt(6, 0.0);

        // Allocate storage for the pressure components (two of these)
        Vector<double> interpolated_p(2, 0.0);

        // Allocate storage for the derivatives of the velocity components
        // w.r.t. local coordinates (s_1 and s_2)
        DenseMatrix<double> interpolated_duds(6, 2, 0.0);

        // Calculate pressure at the integration point
        // -------------------------------------------

        // Loop over pressure degrees of freedom (associated with a single
        // pressure component) in the element
        for (unsigned l = 0; l < n_pres; l++)
        {
          // Cache the shape function
          const double psip_ = psip(l);

          // Loop over the two pressure components
          for (unsigned i = 0; i < 2; i++)
          {
            // Get the value
            const double p_value = this->p_lin_axi_nst(l, i);

            // Add contribution
            interpolated_p[i] += p_value * psip_;
          }
        } // End of loop over the pressure degrees of freedom in the element

        // Calculate eulerian positions, perturbations to these positions,
        // ---------------------------------------------------------------
        // velocities and their derivatives at the integration point
        // ---------------------------------------------------------

        // Loop over the element's nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Cache the shape function
          const double psif_ = psif(l);

          // Loop over the two coordinate directions
          for (unsigned i = 0; i < 2; i++)
          {
            // Calculate the unperturbed position xbar
            interpolated_x[i] += this->raw_nodal_position(l, i) * psif_;

            // Loop over the two coordinate directions (for derivatives)
            for (unsigned j = 0; j < 2; j++)
            {
              interpolated_dxbar_ds(i, j) +=
                this->raw_nodal_position(l, i) * dpsifds(l, j);
            }
          }

          // Loop over the four position perturbations R_k^C, R_k^S, Z_k^C,
          // Z_k^S
          for (unsigned i = 0; i < 4; i++)
          {
            // Get the value
            const double xhat_value =
              this->raw_nodal_value(l, this->xhat_index_lin_axi_nst(l, i));

            // Add contribution
            interpolated_xhat[i] += xhat_value * psif_;

            // Loop over the two coordinate directions (for derivatives)
            for (unsigned j = 0; j < 2; j++)
            {
              interpolated_dxhat_ds(i, j) += xhat_value * dpsifds(l, j);
            }
          }

          // Loop over the six velocity components
          for (unsigned i = 0; i < 6; i++)
          {
            // Get the value
            const double u_value = this->raw_nodal_value(l, u_nodal_index[i]);

            // Add contribution
            interpolated_u[i] += u_value * psif_;

            // Add contribution to dudt
            dudt[i] += this->du_dt_lin_axi_nst(l, i) * psif_;

            // Loop over the two coordinate directions (for derivatives)
            for (unsigned j = 0; j < 2; j++)
            {
              interpolated_duds(i, j) += u_value * dpsifds(l, j);
            }
          }
        } // End of loop over the element's nodes

        // Get the mesh velocity if ALE is enabled
        if (!this->ALE_is_disabled)
        {
          // Loop over the element's nodes
          for (unsigned l = 0; l < n_node; l++)
          {
            // Loop over the two coordinate directions
            for (unsigned i = 0; i < 2; i++)
            {
              // Calculate the derivative of the unperturbed position xbar
              // w.r.t. time
              mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif(l);
            }

            // Loop over the four perturbed positions R_k^C, R_k^S, Z_k^C, Z_k^S
            for (unsigned i = 0; i < 4; i++)
            {
              // Calculate the derivative of the perturbed position xhat
              // w.r.t. time
              dxhat_dt[i] +=
                this->dnodal_position_perturbation_dt_lin_axi_nst(l, i) *
                psif(l);
            }
          }
        }

        // Compute the cosine part of the perturbed jacobian at this ipt
        const double JhatC =
          interpolated_dxbar_ds(0, 0) * interpolated_dxhat_ds(2, 1) +
          interpolated_dxbar_ds(1, 1) * interpolated_dxhat_ds(0, 0) -
          interpolated_dxbar_ds(0, 1) * interpolated_dxhat_ds(2, 0) -
          interpolated_dxbar_ds(1, 0) * interpolated_dxhat_ds(0, 1);

        // Compute the sine part of the perturbed jacobian at this ipt
        const double JhatS =
          interpolated_dxbar_ds(0, 0) * interpolated_dxhat_ds(3, 1) +
          interpolated_dxbar_ds(1, 1) * interpolated_dxhat_ds(1, 0) -
          interpolated_dxbar_ds(0, 1) * interpolated_dxhat_ds(3, 0) -
          interpolated_dxbar_ds(1, 0) * interpolated_dxhat_ds(1, 1);

        // Determine the inverse of the unperturbed jacobian
        const double invJbar = 1.0 / Jbar;

        // Get the user-defined (base flow) body force terms
        Vector<double> body_force(3);
        this->get_body_force_base_flow(time, ipt, interpolated_x, body_force);

        // Get the user-defined (base flow) source function
        const double source =
          this->get_source_base_flow(time, ipt, interpolated_x);

        // Get velocities and their derivatives from base flow problem
        // -----------------------------------------------------------

        // Allocate storage for the velocity components of the base state
        // solution (initialise to zero)
        Vector<double> base_flow_u(3, 0.0);

        // Get the base state solution velocity components
        this->get_base_flow_u(time, ipt, interpolated_x, base_flow_u);

        // Allocate storage for the derivatives of the base state solution's
        // velocity components w.r.t. local coordinates (s_1 and s_2) and
        // global coordinates (r and z)
        // N.B. the derivatives of the base flow components w.r.t. the
        // azimuthal coordinate direction (theta) are always zero since the
        // base flow is axisymmetric
        DenseMatrix<double> base_flow_duds(3, 2, 0.0);
        DenseMatrix<double> base_flow_dudx(3, 2, 0.0);

        // Get the derivatives of the base state solution
        // velocity components w.r.t. local and global coordinates
        this->get_base_flow_duds(time, ipt, interpolated_x, base_flow_duds);
        this->get_base_flow_dudx(time, ipt, interpolated_x, base_flow_dudx);

        // Allocate storage for the base state pressure at the current
        // integration point
        double base_flow_p = 0.0;

        // Allocate storage for the derivatives of the base state solution
        // velocity components w.r.t. time
        Vector<double> base_flow_dudt(3, 0.0);

        // If ALE is enabled, get the base state pressure and the derivatives
        // of the base state velocity w.r.t. time (only needed in this case)
        // if (!ALE_is_disabled)
        {
          this->get_base_flow_p(time, ipt, interpolated_x, base_flow_p);
          this->get_base_flow_dudt(time, ipt, interpolated_x, base_flow_dudt);
        }

        // Compute the following quantities
        const double interpolated_dUdRC =
          invJbar * (interpolated_dxbar_ds(1, 1) * interpolated_duds(0, 0) +
                     base_flow_duds(0, 0) * interpolated_dxhat_ds(2, 1) -
                     interpolated_dxbar_ds(1, 0) * interpolated_duds(0, 1) -
                     base_flow_duds(0, 1) * interpolated_dxhat_ds(2, 0));
        const double interpolated_dUdRS =
          invJbar * (interpolated_dxbar_ds(1, 1) * interpolated_duds(1, 0) +
                     base_flow_duds(0, 0) * interpolated_dxhat_ds(3, 1) -
                     interpolated_dxbar_ds(1, 0) * interpolated_duds(1, 1) -
                     base_flow_duds(0, 1) * interpolated_dxhat_ds(3, 0));
        const double interpolated_dWdRC =
          invJbar * (interpolated_dxbar_ds(1, 1) * interpolated_duds(2, 0) +
                     base_flow_duds(1, 0) * interpolated_dxhat_ds(2, 1) -
                     interpolated_dxbar_ds(1, 0) * interpolated_duds(2, 1) -
                     base_flow_duds(1, 1) * interpolated_dxhat_ds(2, 0));
        const double interpolated_dWdRS =
          invJbar * (interpolated_dxbar_ds(1, 1) * interpolated_duds(3, 0) +
                     base_flow_duds(1, 0) * interpolated_dxhat_ds(3, 1) -
                     interpolated_dxbar_ds(1, 0) * interpolated_duds(3, 1) -
                     base_flow_duds(1, 1) * interpolated_dxhat_ds(3, 0));
        const double interpolated_dVdRC =
          invJbar * (interpolated_dxbar_ds(1, 1) * interpolated_duds(4, 0) +
                     base_flow_duds(2, 0) * interpolated_dxhat_ds(2, 1) -
                     interpolated_dxbar_ds(1, 0) * interpolated_duds(4, 1) -
                     base_flow_duds(2, 1) * interpolated_dxhat_ds(2, 0));
        const double interpolated_dVdRS =
          invJbar * (interpolated_dxbar_ds(1, 1) * interpolated_duds(5, 0) +
                     base_flow_duds(2, 0) * interpolated_dxhat_ds(3, 1) -
                     interpolated_dxbar_ds(1, 0) * interpolated_duds(5, 1) -
                     base_flow_duds(2, 1) * interpolated_dxhat_ds(3, 0));
        const double interpolated_dUdZC =
          invJbar * (interpolated_dxbar_ds(0, 0) * interpolated_duds(0, 1) +
                     base_flow_duds(0, 1) * interpolated_dxhat_ds(0, 0) -
                     interpolated_dxbar_ds(0, 1) * interpolated_duds(0, 0) -
                     base_flow_duds(0, 0) * interpolated_dxhat_ds(0, 1));
        const double interpolated_dUdZS =
          invJbar * (interpolated_dxbar_ds(0, 0) * interpolated_duds(1, 1) +
                     base_flow_duds(0, 1) * interpolated_dxhat_ds(1, 0) -
                     interpolated_dxbar_ds(0, 1) * interpolated_duds(1, 0) -
                     base_flow_duds(0, 0) * interpolated_dxhat_ds(1, 1));
        const double interpolated_dWdZC =
          invJbar * (interpolated_dxbar_ds(0, 0) * interpolated_duds(2, 1) +
                     base_flow_duds(1, 1) * interpolated_dxhat_ds(0, 0) -
                     interpolated_dxbar_ds(0, 1) * interpolated_duds(2, 0) -
                     base_flow_duds(1, 0) * interpolated_dxhat_ds(0, 1));
        const double interpolated_dWdZS =
          invJbar * (interpolated_dxbar_ds(0, 0) * interpolated_duds(3, 1) +
                     base_flow_duds(1, 1) * interpolated_dxhat_ds(1, 0) -
                     interpolated_dxbar_ds(0, 1) * interpolated_duds(3, 0) -
                     base_flow_duds(1, 0) * interpolated_dxhat_ds(1, 1));
        const double interpolated_dVdZC =
          invJbar * (interpolated_dxbar_ds(0, 0) * interpolated_duds(4, 1) +
                     base_flow_duds(2, 1) * interpolated_dxhat_ds(0, 0) -
                     interpolated_dxbar_ds(0, 1) * interpolated_duds(4, 0) -
                     base_flow_duds(2, 0) * interpolated_dxhat_ds(0, 1));
        const double interpolated_dVdZS =
          invJbar * (interpolated_dxbar_ds(0, 0) * interpolated_duds(5, 1) +
                     base_flow_duds(2, 1) * interpolated_dxhat_ds(1, 0) -
                     interpolated_dxbar_ds(0, 1) * interpolated_duds(5, 0) -
                     base_flow_duds(2, 0) * interpolated_dxhat_ds(1, 1));

        // Define the following useful quantities...
        Vector<double> group_A(n_node, 0.0);
        Vector<double> group_B(n_node, 0.0);
        Vector<double> group_C(n_node, 0.0);
        Vector<double> group_D(n_node, 0.0);
        Vector<double> group_E(n_node, 0.0);
        DenseMatrix<double> group_F(n_node, n_node, 0.0);

        // Loop over the element's nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          group_A[l] = interpolated_dxbar_ds(0, 0) * dpsifds(l, 1) -
                       interpolated_dxbar_ds(0, 1) * dpsifds(l, 0);
          group_B[l] = interpolated_dxbar_ds(1, 1) * dpsifds(l, 0) -
                       interpolated_dxbar_ds(1, 0) * dpsifds(l, 1);
          group_C[l] = base_flow_duds(0, 0) * dpsifds(l, 1) -
                       base_flow_duds(0, 1) * dpsifds(l, 0);
          group_D[l] = base_flow_duds(1, 0) * dpsifds(l, 1) -
                       base_flow_duds(1, 1) * dpsifds(l, 0);
          group_E[l] = base_flow_duds(2, 0) * dpsifds(l, 1) -
                       base_flow_duds(2, 1) * dpsifds(l, 0);

          // Loop over the element's nodes again
          for (unsigned l2 = 0; l2 < n_node; l2++)
          {
            group_F(l, l2) =
              dtestfds(l, 0) * dpsifds(l2, 1) - dtestfds(l, 1) * dpsifds(l2, 0);
          }
        }

        // Cache base flow velocities and their derivatives
        const double base_flow_ur = base_flow_u[0];
        const double base_flow_uz = base_flow_u[1];
        const double base_flow_utheta = base_flow_u[2];
        const double base_flow_durdr = base_flow_dudx(0, 0);
        const double base_flow_durdz = base_flow_dudx(0, 1);
        const double base_flow_duzdr = base_flow_dudx(1, 0);
        const double base_flow_duzdz = base_flow_dudx(1, 1);
        const double base_flow_duthetadr = base_flow_dudx(2, 0);
        const double base_flow_duthetadz = base_flow_dudx(2, 1);
        const double base_flow_durdt = base_flow_dudt[0];
        const double base_flow_duzdt = base_flow_dudt[1];
        const double base_flow_duthetadt = base_flow_dudt[2];

        // Cache r-component of position
        const double r = interpolated_x[0];

        // Cache perturbations to nodal positions
        const double interpolated_RC = interpolated_xhat[0];
        const double interpolated_RS = interpolated_xhat[1];
        const double interpolated_ZC = interpolated_xhat[2];
        const double interpolated_ZS = interpolated_xhat[3];

        // Cache temporal derivatives of the perturbations to nodal positions
        const double dRCdt = dxhat_dt[0];
        const double dRSdt = dxhat_dt[1];
        const double dZCdt = dxhat_dt[2];
        const double dZSdt = dxhat_dt[3];

        // Cache unknowns
        const double interpolated_UC = interpolated_u[0];
        const double interpolated_US = interpolated_u[1];
        const double interpolated_WC = interpolated_u[2];
        const double interpolated_WS = interpolated_u[3];
        const double interpolated_VC = interpolated_u[4];
        const double interpolated_VS = interpolated_u[5];
        const double interpolated_PC = interpolated_p[0];
        const double interpolated_PS = interpolated_p[1];

        // Cache temporal derivatives of the unknowns
        const double dUCdt = dudt[0];
        const double dUSdt = dudt[1];
        const double dWCdt = dudt[2];
        const double dWSdt = dudt[3];
        const double dVCdt = dudt[4];
        const double dVSdt = dudt[5];

        // ==================
        // MOMENTUM EQUATIONS
        // ==================

        // Loop over the fluid test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Cache test functions and their derivatives
          const double testf_ = testf(l);
          const double dtestfdr = dtestfdx(l, 0);
          const double dtestfdz = dtestfdx(l, 1);

          // Compute the following useful quantities...
          const double dtestfdRC =
            invJbar * (dtestfds(l, 0) * interpolated_dxhat_ds(2, 1) -
                       dtestfds(l, 1) * interpolated_dxhat_ds(2, 0));
          const double dtestfdRS =
            invJbar * (dtestfds(l, 0) * interpolated_dxhat_ds(3, 1) -
                       dtestfds(l, 1) * interpolated_dxhat_ds(3, 0));
          const double dtestfdZC =
            invJbar * (dtestfds(l, 1) * interpolated_dxhat_ds(0, 0) -
                       dtestfds(l, 0) * interpolated_dxhat_ds(0, 1));
          const double dtestfdZS =
            invJbar * (dtestfds(l, 1) * interpolated_dxhat_ds(1, 0) -
                       dtestfds(l, 0) * interpolated_dxhat_ds(1, 1));

          // ---------------------------------------------
          // FIRST (RADIAL) MOMENTUM EQUATION: COSINE PART
          // ---------------------------------------------

          // Get local equation number of first velocity value at this node
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[0]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -=
              scaled_re_st * r * dUCdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re_st * interpolated_RC *
                                    base_flow_durdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_ur *
                                    interpolated_dUdRC * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[0] *
                                    interpolated_dUdRC * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_UC *
                                    base_flow_durdr * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dRCdt * base_flow_durdr * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RC * base_flow_ur *
                                    base_flow_durdr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RC *
                                    mesh_velocity[0] * base_flow_durdr *
                                    testf_ * Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    interpolated_US * testf_ * Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    base_flow_durdr * interpolated_RS * testf_ *
                                    Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    base_flow_durdz * interpolated_ZS * testf_ *
                                    Jbar * w;
            residuals[local_eqn] += 2 * scaled_re * base_flow_utheta *
                                    interpolated_VC * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_uz *
                                    interpolated_dUdZC * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[1] *
                                    interpolated_dUdZC * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_WC *
                                    base_flow_durdz * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dZCdt * base_flow_durdz * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RC * base_flow_uz *
                                    base_flow_durdz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RC *
                                    mesh_velocity[1] * base_flow_durdz *
                                    testf_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RC * body_force[0] * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * interpolated_RC * G[0] * testf_ * Jbar * w;
            residuals[local_eqn] += r * base_flow_p * dtestfdRC * Jbar * w;
            residuals[local_eqn] += r * interpolated_PC * dtestfdr * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RC * base_flow_p * dtestfdr * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * r *
                                    base_flow_durdr * dtestfdRC * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * r *
                                    interpolated_dUdRC * dtestfdr * Jbar * w;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * r *
                                    base_flow_durdr * JhatC * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) *
                                    interpolated_RC * base_flow_durdr *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    base_flow_duthetadr * dtestfdr *
                                    interpolated_RS * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    base_flow_duthetadr * dtestfdz *
                                    interpolated_ZS * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    interpolated_dVdRS * testf_ * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * k * k * interpolated_UC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_durdr *
                                    interpolated_RC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_durdz *
                                    interpolated_ZC * testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * k * base_flow_utheta *
                                    dtestfdr * interpolated_RS * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * k * base_flow_utheta *
                                    dtestfdz * interpolated_ZS * Jbar * w / r;
            residuals[local_eqn] -=
              visc_ratio * k * interpolated_VS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * interpolated_RS *
                                    base_flow_utheta * testf_ * Jbar * w /
                                    (r * r);
            residuals[local_eqn] -= visc_ratio * this->Gamma[0] * r *
                                    base_flow_duzdr * dtestfdZC * Jbar * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[0] * r *
                                    interpolated_dWdRC * dtestfdz * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] * r *
                                    base_flow_duzdr * JhatC * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[0] *
                                    interpolated_RC * base_flow_duzdr *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_durdz * dtestfdZC * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dUdZC * dtestfdz * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_durdz * JhatC * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RC *
                                    base_flow_durdz * dtestfdz * Jbar * w;
            residuals[local_eqn] += interpolated_PC * testf_ * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    interpolated_VS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_duthetadr * interpolated_RS *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_duthetadz * interpolated_ZS *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) *
                                    interpolated_UC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) *
                                    base_flow_ur * interpolated_RC * testf_ *
                                    Jbar * w / (r * r);
            residuals[local_eqn] -=
              scaled_re_st * r * base_flow_durdt * testf_ * JhatC * w;
            residuals[local_eqn] += scaled_re * base_flow_utheta *
                                    base_flow_utheta * testf_ * JhatC * w;
            residuals[local_eqn] += r * body_force[0] * testf_ * JhatC * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * r * G[0] * testf_ * JhatC * w;
            residuals[local_eqn] -=
              k * visc_ratio * base_flow_utheta * testf_ * JhatS * w / r;
            residuals[local_eqn] += base_flow_p * testf_ * JhatC * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) *
                                    base_flow_ur * testf_ * JhatC * w / r;

            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[0]);
                if (local_unknown >= 0)
                {
                  if (flag == 2)
                  {
                    // Add the mass matrix
                    mass_matrix(local_eqn, local_unknown) +=
                      scaled_re_st * r * psif[l2] * testf_ * Jbar * w;
                  }

                  // Add contributions to the Jacobian matrix
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_B[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_durdr * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_uz * group_A[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[1] * group_A[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * group_B[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * k * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_A[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * psif[l2] * testf_ *
                    Jbar * w / r;
                }

                // Radial velocity component (sine part) U_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Axial velocity component (cosine part) W_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_durdz * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[0] * r * group_B[l2] * dtestfdz *
                    w;
                }

                // Axial velocity component (sine part) W_k^S
                // has no contribution

                // Azimuthal velocity component (cosine part) V_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[4]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    2.0 * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Azimuthal velocity component (sine part) V_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[5]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * psif[l2] *
                    testf_ * Jbar * w / r;
                }

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * psif[l2] * base_flow_durdt * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_durdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * base_flow_durdr *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[0] *
                    base_flow_durdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * r * base_flow_uz * group_C[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * mesh_velocity[1] * group_C[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_uz * base_flow_durdz *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[1] *
                    base_flow_durdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * body_force[0] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * psif[l2] * G[0] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_p * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * base_flow_durdr *
                    group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * psif[l2] *
                    base_flow_durdr * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_durdr * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * r * base_flow_duzdr *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * r * base_flow_duzdr *
                    group_B[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[0] * psif[l2] * base_flow_duzdr *
                    dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_durdz * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * group_C[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_durdz * group_B[l2] * dtestfdz *
                    w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_durdz *
                                                        dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    psif[l2] * testf_ * Jbar * w / (r * r);
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_durdt * testf_ * group_B[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * base_flow_utheta * base_flow_utheta * testf_ *
                    group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[0] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[0] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_p * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_B[l2] * w / r;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * base_flow_durdr *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * base_flow_duthetadr *
                    dtestfdr * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * base_flow_utheta * dtestfdr * psif[l2] *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * psif[l2] * base_flow_utheta * testf_ *
                    Jbar * w / (r * r);
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k *
                    base_flow_duthetadr * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * base_flow_utheta * testf_ * group_B[l2] *
                    w / r;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_C[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_C[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_durdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * base_flow_p * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * base_flow_durdr *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * group_C[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * base_flow_durdr *
                    group_A[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_durdz * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[0] * r * group_D[l2] * dtestfdz *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * r * base_flow_duzdr *
                    group_A[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_durdz * group_A[l2] * dtestfdz *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_durdt * testf_ * group_A[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * base_flow_utheta * base_flow_utheta * testf_ *
                    group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[0] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[0] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_p * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_A[l2] * w / r;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * base_flow_durdz *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * base_flow_duthetadr *
                    dtestfdz * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * base_flow_utheta * dtestfdz * psif[l2] *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k *
                    base_flow_duthetadz * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * base_flow_utheta * testf_ * group_A[l2] *
                    w / r;
                }

              } // End of loop over velocity shape functions

              // Now loop over pressure shape functions
              // (This is the contribution from pressure gradient)
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                // Cosine part P_k^C
                local_unknown = this->p_local_eqn(l2, 0);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * psip[l2] * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psip[l2] * testf_ * Jbar * w;
                }

                // Sine part P_k^S has no contribution

              } // End of loop over pressure shape functions

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

          // --------------------------------------------
          // SECOND (RADIAL) MOMENTUM EQUATION: SINE PART
          // --------------------------------------------

          // Get local equation number of second velocity value at this node
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[1]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -=
              scaled_re_st * r * dUSdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re_st * interpolated_RS *
                                    base_flow_durdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_ur *
                                    interpolated_dUdRS * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[0] *
                                    interpolated_dUdRS * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_US *
                                    base_flow_durdr * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dRSdt * base_flow_durdr * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RS * base_flow_ur *
                                    base_flow_durdr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RS *
                                    mesh_velocity[0] * base_flow_durdr *
                                    testf_ * Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    interpolated_UC * testf_ * Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    base_flow_durdr * interpolated_RC * testf_ *
                                    Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    base_flow_durdz * interpolated_ZC * testf_ *
                                    Jbar * w;
            residuals[local_eqn] += 2 * scaled_re * base_flow_utheta *
                                    interpolated_VS * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_uz *
                                    interpolated_dUdZS * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[1] *
                                    interpolated_dUdZS * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_WS *
                                    base_flow_durdz * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dZSdt * base_flow_durdz * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RS * base_flow_uz *
                                    base_flow_durdz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RS *
                                    mesh_velocity[1] * base_flow_durdz *
                                    testf_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RS * body_force[0] * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * interpolated_RS * G[0] * testf_ * Jbar * w;
            residuals[local_eqn] += r * base_flow_p * dtestfdRS * Jbar * w;
            residuals[local_eqn] += r * interpolated_PS * dtestfdr * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RS * base_flow_p * dtestfdr * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * r *
                                    base_flow_durdr * dtestfdRS * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * r *
                                    interpolated_dUdRS * dtestfdr * Jbar * w;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * r *
                                    base_flow_durdr * JhatS * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) *
                                    interpolated_RS * base_flow_durdr *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    base_flow_duthetadr * dtestfdr *
                                    interpolated_RC * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    base_flow_duthetadr * dtestfdz *
                                    interpolated_ZC * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    interpolated_dVdRC * testf_ * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * k * k * interpolated_US * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_durdr *
                                    interpolated_RS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_durdz *
                                    interpolated_ZS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * base_flow_utheta *
                                    dtestfdr * interpolated_RC * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * base_flow_utheta *
                                    dtestfdz * interpolated_ZC * Jbar * w / r;
            residuals[local_eqn] +=
              visc_ratio * k * interpolated_VC * testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * k * interpolated_RC *
                                    base_flow_utheta * testf_ * Jbar * w /
                                    (r * r);
            residuals[local_eqn] -= visc_ratio * this->Gamma[0] * r *
                                    base_flow_duzdr * dtestfdZS * Jbar * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[0] * r *
                                    interpolated_dWdRS * dtestfdz * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] * r *
                                    base_flow_duzdr * JhatS * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[0] *
                                    interpolated_RS * base_flow_duzdr *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_durdz * dtestfdZS * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dUdZS * dtestfdz * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_durdz * JhatS * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RS *
                                    base_flow_durdz * dtestfdz * Jbar * w;
            residuals[local_eqn] += interpolated_PS * testf_ * Jbar * w;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    interpolated_VC * testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_duthetadr * interpolated_RC *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_duthetadz * interpolated_ZC *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) *
                                    interpolated_US * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) *
                                    base_flow_ur * interpolated_RS * testf_ *
                                    Jbar * w / (r * r);
            residuals[local_eqn] -=
              scaled_re_st * r * base_flow_durdt * testf_ * JhatS * w;
            residuals[local_eqn] += scaled_re * base_flow_utheta *
                                    base_flow_utheta * testf_ * JhatS * w;
            residuals[local_eqn] += r * body_force[0] * testf_ * JhatS * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * r * G[0] * testf_ * JhatS * w;
            residuals[local_eqn] +=
              k * visc_ratio * base_flow_utheta * testf_ * JhatC * w / r;
            residuals[local_eqn] += base_flow_p * testf_ * JhatS * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) *
                                    base_flow_ur * testf_ * JhatS * w / r;


            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[0]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Radial velocity component (sine part) U_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  if (flag == 2)
                  {
                    // Add the mass matrix
                    mass_matrix(local_eqn, local_unknown) +=
                      scaled_re_st * r * psif[l2] * testf_ * Jbar * w;
                  }

                  // Add contributions to the Jacobian matrix
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_B[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_durdr * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_uz * group_A[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[1] * group_A[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * group_B[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * k * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_A[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * psif[l2] * testf_ *
                    Jbar * w / r;
                }

                // Axial velocity component (cosine part) W_k^C
                // has no contribution

                // Axial velocity component (sine part) W_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[3]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_durdz * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[0] * r * group_B[l2] * dtestfdz *
                    w;
                }

                // Azimuthal velocity component (cosine part) V_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[4]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * psif[l2] *
                    testf_ * Jbar * w / r;
                }

                // Azimuthal velocity component (sine part) V_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[5]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    2.0 * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * base_flow_durdr *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * base_flow_duthetadr *
                    dtestfdr * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * base_flow_utheta * dtestfdr * psif[l2] *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * psif[l2] * base_flow_utheta * testf_ *
                    Jbar * w / (r * r);
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k *
                    base_flow_duthetadr * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * base_flow_utheta * testf_ * group_B[l2] *
                    w / r;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * psif[l2] * base_flow_durdt * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_durdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * base_flow_durdr *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[0] *
                    base_flow_durdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * r * base_flow_uz * group_C[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * mesh_velocity[1] * group_C[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_uz * base_flow_durdz *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[1] *
                    base_flow_durdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * body_force[0] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * psif[l2] * G[0] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_p * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * base_flow_durdr *
                    group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * psif[l2] *
                    base_flow_durdr * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_durdr * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * r * base_flow_duzdr *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * r * base_flow_duzdr *
                    group_B[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[0] * psif[l2] * base_flow_duzdr *
                    dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_durdz * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * group_C[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_durdz * group_B[l2] * dtestfdz *
                    w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_durdz *
                                                        dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    psif[l2] * testf_ * Jbar * w / (r * r);
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_durdt * testf_ * group_B[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * base_flow_utheta * base_flow_utheta * testf_ *
                    group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[0] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[0] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_p * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_B[l2] * w / r;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * base_flow_durdz *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * base_flow_duthetadr *
                    dtestfdz * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * base_flow_utheta * dtestfdz * psif[l2] *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k *
                    base_flow_duthetadz * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * base_flow_utheta * testf_ * group_A[l2] *
                    w / r;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_C[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_C[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_durdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * base_flow_p * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * base_flow_durdr *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * group_C[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * r * base_flow_durdr *
                    group_A[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_durdz * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * r * base_flow_duzdr *
                    dtestfdz * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[0] * r * group_D[l2] * dtestfdz *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_durdz * group_A[l2] * dtestfdz *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_durdt * testf_ * group_A[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * base_flow_utheta * base_flow_utheta * testf_ *
                    group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[0] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[0] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_p * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_A[l2] * w / r;
                }

              } // End of loop over velocity shape functions

              // Now loop over pressure shape functions
              // (This is the contribution from pressure gradient)
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                // Cosine part P_k^C has no contribution

                // Sine part P_k^S
                local_unknown = this->p_local_eqn(l2, 1);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * psip[l2] * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psip[l2] * testf_ * Jbar * w;
                }
              } // End of loop over pressure shape functions

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

          // --------------------------------------------
          // THIRD (AXIAL) MOMENTUM EQUATION: COSINE PART
          // --------------------------------------------

          // Get local equation number of third velocity value at this node
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[2]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -=
              scaled_re_st * r * dWCdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re_st * interpolated_RC *
                                    base_flow_duzdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_ur *
                                    interpolated_dWdRC * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[0] *
                                    interpolated_dWdRC * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_UC *
                                    base_flow_duzdr * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dRCdt * base_flow_duzdr * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RC * base_flow_ur *
                                    base_flow_duzdr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RC *
                                    mesh_velocity[0] * base_flow_duzdr *
                                    testf_ * Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    interpolated_WS * testf_ * Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    base_flow_duzdr * interpolated_RS * testf_ *
                                    Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    base_flow_duzdz * interpolated_ZS * testf_ *
                                    Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_uz *
                                    interpolated_dWdZC * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[1] *
                                    interpolated_dWdZC * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_WC *
                                    base_flow_duzdz * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dZCdt * base_flow_duzdz * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RC * base_flow_uz *
                                    base_flow_duzdz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RC *
                                    mesh_velocity[1] * base_flow_duzdz *
                                    testf_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RC * body_force[1] * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * interpolated_RC * G[1] * testf_ * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_duzdr * dtestfdRC * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dWdRC * dtestfdr * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_duzdr * JhatC * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RC *
                                    base_flow_duzdr * dtestfdr * Jbar * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[1] * r *
                                    base_flow_durdz * dtestfdRC * Jbar * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[1] * r *
                                    interpolated_dUdZC * dtestfdr * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[1] * r *
                                    base_flow_durdz * JhatC * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[1] *
                                    interpolated_RC * base_flow_durdz *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * k * k * interpolated_WC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_duzdr *
                                    interpolated_RC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_duzdz *
                                    interpolated_ZC * testf_ * Jbar * w / r;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[1] *
                                    base_flow_duthetadz * dtestfdr *
                                    interpolated_RS * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[1] *
                                    base_flow_duthetadz * dtestfdz *
                                    interpolated_ZS * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[1] *
                                    interpolated_dVdZS * testf_ * Jbar * w;
            residuals[local_eqn] += r * base_flow_p * dtestfdZC * Jbar * w;
            residuals[local_eqn] += r * interpolated_PC * dtestfdz * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RC * base_flow_p * dtestfdz * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[1]) * r *
                                    base_flow_duzdz * dtestfdZC * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[1]) * r *
                                    interpolated_dWdZC * dtestfdz * Jbar * w;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[1]) * r *
                                    base_flow_duzdz * JhatC * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[1]) *
                                    interpolated_RC * base_flow_duzdz *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] -=
              scaled_re_st * r * base_flow_duzdt * testf_ * JhatC * w;
            residuals[local_eqn] += r * body_force[1] * testf_ * JhatC * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * r * G[1] * testf_ * JhatC * w;


            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[0]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duzdr * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[1] * r * group_A[l2] * dtestfdr *
                    w;
                }

                // Radial velocity component (sine part) U_k^S
                // has no contribution

                // Axial velocity component (cosine part) W_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  if (flag == 2)
                  {
                    // Add the mass matrix
                    mass_matrix(local_eqn, local_unknown) +=
                      scaled_re_st * r * psif[l2] * testf_ * Jbar * w;
                  }

                  // Add contributions to the Jacobian matrix
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_B[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_uz * group_A[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[1] * group_A[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duzdz * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * k * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * group_A[l2] *
                    dtestfdz * w;
                }

                // Axial velocity component (sine part) W_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[3]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Azimuthal velocity component (cosine part) V_k^C
                // has no contribution

                // Azimuthal velocity component (sine part) V_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[5]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[1] * group_A[l2] * testf_ * w;
                }

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * psif[l2] * base_flow_duzdt * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duzdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * base_flow_duzdr *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[0] *
                    base_flow_duzdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * r * base_flow_uz * group_D[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * mesh_velocity[1] * group_D[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_uz * base_flow_duzdz *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[1] *
                    base_flow_duzdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * body_force[1] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * psif[l2] * G[1] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duzdr * group_B[l2] * dtestfdr *
                    w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_duzdr *
                                                        dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[1] * r * group_C[l2] * dtestfdr *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[1] * r * base_flow_durdz *
                    group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[1] * psif[l2] * base_flow_durdz *
                    dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_duzdr * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    r * base_flow_p * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_p * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * base_flow_duzdz *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * group_D[l2] *
                    dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * base_flow_duzdz *
                    group_B[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[1]) * psif[l2] *
                    base_flow_duzdz * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duzdt * testf_ * group_B[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[1] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[1] * testf_ * group_B[l2] * w;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * base_flow_duzdr *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[1] * base_flow_duthetadz *
                    dtestfdr * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[1] * group_E[l2] * testf_ * w;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_D[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_D[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duzdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * base_flow_duzdr * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_D[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duzdr * group_A[l2] * dtestfdr *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[1] * r * base_flow_durdz *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[1] * r * base_flow_durdz *
                    group_A[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_duzdz * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * base_flow_duzdz *
                    group_A[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duzdt * testf_ * group_A[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[1] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[1] * testf_ * group_A[l2] * w;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * base_flow_duzdz *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[1] * base_flow_duthetadz *
                    dtestfdz * psif[l2] * Jbar * w;
                }

              } // End of loop over velocity shape functions

              // Now loop over pressure shape functions
              // (This is the contribution from pressure gradient)
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                // Cosine part P_k^C
                local_unknown = this->p_local_eqn(l2, 0);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * psip[l2] * dtestfdz * Jbar * w;
                }

                // Sine part P_k^S has no contribution

              } // End of loop over pressure shape functions

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

          // -------------------------------------------
          // FOURTH (AXIAL) MOMENTUM EQUATION: SINE PART
          // -------------------------------------------

          // Get local equation number of fourth velocity value at this node
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[3]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -=
              scaled_re_st * r * dWSdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re_st * interpolated_RS *
                                    base_flow_duzdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_ur *
                                    interpolated_dWdRS * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[0] *
                                    interpolated_dWdRS * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_US *
                                    base_flow_duzdr * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dRSdt * base_flow_duzdr * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RS * base_flow_ur *
                                    base_flow_duzdr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RS *
                                    mesh_velocity[0] * base_flow_duzdr *
                                    testf_ * Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    interpolated_WC * testf_ * Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    base_flow_duzdr * interpolated_RC * testf_ *
                                    Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    base_flow_duzdz * interpolated_ZC * testf_ *
                                    Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_uz *
                                    interpolated_dWdZS * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[1] *
                                    interpolated_dWdZS * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_WS *
                                    base_flow_duzdz * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_st * r * dZSdt * base_flow_duzdz * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RS * base_flow_uz *
                                    base_flow_duzdz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RS *
                                    mesh_velocity[1] * base_flow_duzdz *
                                    testf_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RS * body_force[1] * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * interpolated_RS * G[1] * testf_ * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_duzdr * dtestfdRS * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dWdRS * dtestfdr * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_duzdr * JhatS * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RS *
                                    base_flow_duzdr * dtestfdr * Jbar * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[1] * r *
                                    base_flow_durdz * dtestfdRS * Jbar * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[1] * r *
                                    interpolated_dUdZS * dtestfdr * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[1] * r *
                                    base_flow_durdz * JhatS * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * this->Gamma[1] *
                                    interpolated_RS * base_flow_durdz *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * k * k * interpolated_WS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_duzdr *
                                    interpolated_RS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * k * base_flow_duzdz *
                                    interpolated_ZS * testf_ * Jbar * w / r;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[1] *
                                    base_flow_duthetadz * dtestfdr *
                                    interpolated_RC * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[1] *
                                    base_flow_duthetadz * dtestfdz *
                                    interpolated_ZC * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[1] *
                                    interpolated_dVdZC * testf_ * Jbar * w;
            residuals[local_eqn] += r * base_flow_p * dtestfdZS * Jbar * w;
            residuals[local_eqn] += r * interpolated_PS * dtestfdz * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RS * base_flow_p * dtestfdz * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[1]) * r *
                                    base_flow_duzdz * dtestfdZS * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[1]) * r *
                                    interpolated_dWdZS * dtestfdz * Jbar * w;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[1]) * r *
                                    base_flow_duzdz * JhatS * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[1]) *
                                    interpolated_RS * base_flow_duzdz *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] -=
              scaled_re_st * r * base_flow_duzdt * testf_ * JhatS * w;
            residuals[local_eqn] += r * body_force[1] * testf_ * JhatS * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * r * G[1] * testf_ * JhatS * w;


            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                // has no contribution

                // Radial velocity component (sine part) U_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duzdr * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[1] * r * group_A[l2] * dtestfdr *
                    w;
                }

                // Axial velocity component (cosine part) W_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Axial velocity component (sine part) W_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[3]);
                if (local_unknown >= 0)
                {
                  if (flag == 2)
                  {
                    // Add the mass matrix
                    mass_matrix(local_eqn, local_unknown) +=
                      scaled_re_st * r * psif[l2] * testf_ * Jbar * w;
                  }

                  // Add contributions to the Jacobian matrix
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_B[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_uz * group_A[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[1] * group_A[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duzdz * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * k * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * group_A[l2] *
                    dtestfdz * w;
                }

                // Azimuthal velocity component (cosine part) V_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[4]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[1] * group_A[l2] * testf_ * w;
                }

                // Azimuthal velocity component (sine part) V_k^S
                // has no contribution

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * base_flow_duzdr *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[1] * base_flow_duthetadz *
                    dtestfdr * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[1] * group_E[l2] * testf_ * w;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * psif[l2] * base_flow_duzdt * testf_ * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duzdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * base_flow_duzdr *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[0] *
                    base_flow_duzdr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * r * base_flow_uz * group_D[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * mesh_velocity[1] * group_D[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_uz * base_flow_duzdz *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[1] *
                    base_flow_duzdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * body_force[1] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * psif[l2] * G[1] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duzdr * group_B[l2] * dtestfdr *
                    w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_duzdr *
                                                        dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[1] * r * group_C[l2] * dtestfdr *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[1] * r * base_flow_durdz *
                    group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[1] * psif[l2] * base_flow_durdz *
                    dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_duzdr * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    r * base_flow_p * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_p * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * base_flow_duzdz *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * group_D[l2] *
                    dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * base_flow_duzdz *
                    group_B[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[1]) * psif[l2] *
                    base_flow_duzdz * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duzdt * testf_ * group_B[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[1] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[1] * testf_ * group_B[l2] * w;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * base_flow_duzdz *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[1] * base_flow_duthetadz *
                    dtestfdz * psif[l2] * Jbar * w;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_D[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_D[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duzdz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * base_flow_duzdr * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_D[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duzdr * group_A[l2] * dtestfdr *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * this->Gamma[1] * r * base_flow_durdz *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[1] * r * base_flow_durdz *
                    group_A[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * k * base_flow_duzdz * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[1]) * r * base_flow_duzdz *
                    group_A[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duzdt * testf_ * group_A[l2] *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[1] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[1] * testf_ * group_A[l2] * w;
                }

              } // End of loop over velocity shape functions

              // Now loop over pressure shape functions
              // (This is the contribution from pressure gradient)
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                // Cosine part P_k^C has no contribution

                // Sine part P_k^S
                local_unknown = this->p_local_eqn(l2, 1);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * psip[l2] * dtestfdz * Jbar * w;
                }
              } // End of loop over pressure shape functions

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

          // ------------------------------------------------
          // FIFTH (AZIMUTHAL) MOMENTUM EQUATION: COSINE PART
          // ------------------------------------------------

          // Get local equation number of fifth velocity value at this node
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[4]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -=
              scaled_re_st * r * dVCdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re_st * interpolated_RC *
                                    base_flow_duthetadt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_ur *
                                    interpolated_dVdRC * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[0] *
                                    interpolated_dVdRC * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_UC *
                                    base_flow_duthetadr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * dRCdt *
                                    base_flow_duthetadr * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RC * base_flow_ur *
                                    base_flow_duthetadr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RC *
                                    mesh_velocity[0] * base_flow_duthetadr *
                                    testf_ * Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    interpolated_VS * testf_ * Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    base_flow_duthetadr * interpolated_RS *
                                    testf_ * Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    base_flow_duthetadz * interpolated_ZS *
                                    testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * base_flow_utheta *
                                    interpolated_UC * testf_ * Jbar * w;
            residuals[local_eqn] -=
              scaled_re * interpolated_VC * base_flow_ur * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_uz *
                                    interpolated_dVdZC * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[1] *
                                    interpolated_dVdZC * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_WC *
                                    base_flow_duthetadz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * dZCdt *
                                    base_flow_duthetadz * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RC * base_flow_uz *
                                    base_flow_duthetadz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RC *
                                    mesh_velocity[1] * base_flow_duthetadz *
                                    testf_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RC * body_force[2] * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * interpolated_RC * G[2] * testf_ * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_duthetadr * dtestfdRC * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dVdRC * dtestfdr * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_duthetadr * JhatC * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RC *
                                    base_flow_duthetadr * dtestfdr * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    interpolated_US * dtestfdr * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    base_flow_durdr * interpolated_RS *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    base_flow_durdz * interpolated_ZS *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] *
                                    base_flow_utheta * dtestfdRC * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] *
                                    interpolated_VC * dtestfdr * Jbar * w;
            residuals[local_eqn] -=
              k * base_flow_p * dtestfdr * interpolated_RS * Jbar * w;
            residuals[local_eqn] -=
              k * base_flow_p * dtestfdz * interpolated_ZS * Jbar * w;
            residuals[local_eqn] -= k * interpolated_PS * testf_ * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    k * interpolated_VC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    k * base_flow_duthetadr * interpolated_RC *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    k * base_flow_duthetadz * interpolated_ZC *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_ur * dtestfdr * interpolated_RS *
                                    Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_ur * dtestfdz * interpolated_ZS *
                                    Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    interpolated_US * testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    interpolated_RS * base_flow_ur * testf_ *
                                    Jbar * w / (r * r);
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    interpolated_WS * dtestfdz * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    base_flow_duzdr * interpolated_RS *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    base_flow_duzdz * interpolated_ZS *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_duthetadz * dtestfdZC * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dVdZC * dtestfdz * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_duthetadz * JhatC * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RC *
                                    base_flow_duthetadz * dtestfdz * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] *
                                    interpolated_dVdRC * testf_ * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * k * interpolated_US * testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * k * base_flow_durdr *
                                    interpolated_RS * testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * k * base_flow_durdz *
                                    interpolated_ZS * testf_ * Jbar * w / r;
            residuals[local_eqn] -=
              visc_ratio * interpolated_VC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * base_flow_utheta *
                                    interpolated_RC * testf_ * Jbar * w /
                                    (r * r);
            residuals[local_eqn] -=
              scaled_re_st * r * base_flow_duthetadt * testf_ * JhatC * w;
            residuals[local_eqn] -=
              scaled_re * base_flow_utheta * base_flow_ur * testf_ * JhatC * w;
            residuals[local_eqn] += r * body_force[2] * testf_ * JhatC * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * r * G[2] * testf_ * JhatC * w;
            residuals[local_eqn] -= k * base_flow_p * testf_ * JhatS * w;
            residuals[local_eqn] += k * visc_ratio * (1.0 + this->Gamma[0]) *
                                    base_flow_ur * testf_ * JhatS * w / r;
            residuals[local_eqn] -=
              visc_ratio * base_flow_utheta * testf_ * JhatC * w / r;


            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[0]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duthetadr * testf_ *
                    Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * base_flow_utheta * psif[l2] * testf_ * Jbar * w;
                }

                // Radial velocity component (sine part) U_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * psif[l2] * dtestfdr *
                    Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * psif[l2] *
                    testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * psif[l2] * testf_ * Jbar * w / r;
                }

                // Axial velocity component (cosine part) W_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duthetadz * testf_ *
                    Jbar * w;
                }

                // Axial velocity component (sine part) W_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[3]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * psif[l2] * dtestfdz *
                    Jbar * w;
                }

                // Azimuthal velocity component (cosine part) V_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[4]);
                if (local_unknown >= 0)
                {
                  if (flag == 2)
                  {
                    // Add the mass matrix
                    mass_matrix(local_eqn, local_unknown) +=
                      scaled_re_st * r * psif[l2] * testf_ * Jbar * w;
                  }

                  // Add contributions to the Jacobian matrix
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_B[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_uz * group_A[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[1] * group_A[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * psif[l2] * dtestfdr * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * k * psif[l2] *
                    testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_A[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * psif[l2] * testf_ * Jbar * w / r;
                }

                // Azimuthal velocity component (sine part) V_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[5]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * psif[l2] * base_flow_duthetadt * testf_ *
                    Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duthetadr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * base_flow_duthetadr *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[0] *
                    base_flow_duthetadr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * r * base_flow_uz * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * mesh_velocity[1] * group_E[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_uz * base_flow_duthetadz *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[1] *
                    base_flow_duthetadz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * body_force[2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * psif[l2] * G[2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadr * group_B[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_duthetadr *
                                                        dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * k *
                    base_flow_duthetadr * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadz * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * group_E[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadz * group_B[l2] *
                    dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_duthetadz *
                                                        dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * base_flow_utheta * psif[l2] * testf_ * Jbar *
                    w / (r * r);
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duthetadt * testf_ *
                    group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * base_flow_utheta * base_flow_ur * testf_ *
                    group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[2] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[2] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * base_flow_utheta * testf_ * group_B[l2] * w /
                    r;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * base_flow_duthetadr *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * base_flow_durdr *
                    psif[l2] * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * base_flow_p * dtestfdr * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * base_flow_ur *
                    dtestfdr * psif[l2] * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * psif[l2] *
                    base_flow_ur * testf_ * Jbar * w / (r * r);
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * base_flow_duzdr *
                    psif[l2] * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * base_flow_durdr * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    k * base_flow_p * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_B[l2] * w / r;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_E[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duthetadz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * base_flow_duthetadr * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_E[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadr * group_A[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * base_flow_utheta *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * k *
                    base_flow_duthetadz * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadz * group_A[l2] *
                    dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duthetadt * testf_ *
                    group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * base_flow_utheta * base_flow_ur * testf_ *
                    group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[2] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[2] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * base_flow_utheta * testf_ * group_A[l2] * w /
                    r;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * base_flow_duthetadz *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * base_flow_durdz *
                    psif[l2] * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * base_flow_p * dtestfdz * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * base_flow_ur *
                    dtestfdz * psif[l2] * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * base_flow_duzdz *
                    psif[l2] * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * base_flow_durdz * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    k * base_flow_p * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_A[l2] * w / r;
                }

              } // End of loop over velocity shape functions

              // Now loop over pressure shape functions
              // (This is the contribution from pressure gradient)
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                // Cosine part P_k^C has no contribution

                // Sine part P_k^S
                local_unknown = this->p_local_eqn(l2, 1);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * psip[l2] * testf_ * Jbar * w;
                }
              } // End of loop over pressure shape functions

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

          // ----------------------------------------------
          // SIXTH (AZIMUTHAL) MOMENTUM EQUATION: SINE PART
          // ----------------------------------------------

          // Get local equation number of sixth velocity value at this node
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[5]);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -=
              scaled_re_st * r * dVSdt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re_st * interpolated_RS *
                                    base_flow_duthetadt * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_ur *
                                    interpolated_dVdRS * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[0] *
                                    interpolated_dVdRS * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_US *
                                    base_flow_duthetadr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * dRSdt *
                                    base_flow_duthetadr * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RS * base_flow_ur *
                                    base_flow_duthetadr * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RS *
                                    mesh_velocity[0] * base_flow_duthetadr *
                                    testf_ * Jbar * w;
            residuals[local_eqn] += k * scaled_re * base_flow_utheta *
                                    interpolated_VC * testf_ * Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    base_flow_duthetadr * interpolated_RC *
                                    testf_ * Jbar * w;
            residuals[local_eqn] -= k * scaled_re * base_flow_utheta *
                                    base_flow_duthetadz * interpolated_ZC *
                                    testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * base_flow_utheta *
                                    interpolated_US * testf_ * Jbar * w;
            residuals[local_eqn] -=
              scaled_re * interpolated_VS * base_flow_ur * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * base_flow_uz *
                                    interpolated_dVdZS * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * mesh_velocity[1] *
                                    interpolated_dVdZS * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * r * interpolated_WS *
                                    base_flow_duthetadz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * r * dZSdt *
                                    base_flow_duthetadz * testf_ * Jbar * w;
            residuals[local_eqn] -= scaled_re * interpolated_RS * base_flow_uz *
                                    base_flow_duthetadz * testf_ * Jbar * w;
            residuals[local_eqn] += scaled_re_st * interpolated_RS *
                                    mesh_velocity[1] * base_flow_duthetadz *
                                    testf_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RS * body_force[2] * testf_ * Jbar * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * interpolated_RS * G[2] * testf_ * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_duthetadr * dtestfdRS * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dVdRS * dtestfdr * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_duthetadr * JhatS * dtestfdr * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RS *
                                    base_flow_duthetadr * dtestfdr * Jbar * w;
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    interpolated_UC * dtestfdr * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    base_flow_durdr * interpolated_RC *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    base_flow_durdz * interpolated_ZC *
                                    dtestfdr * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] *
                                    base_flow_utheta * dtestfdRS * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] *
                                    interpolated_VS * dtestfdr * Jbar * w;
            residuals[local_eqn] +=
              k * base_flow_p * dtestfdr * interpolated_RC * Jbar * w;
            residuals[local_eqn] +=
              k * base_flow_p * dtestfdz * interpolated_ZC * Jbar * w;
            residuals[local_eqn] += k * interpolated_PC * testf_ * Jbar * w;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    k * interpolated_VS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    k * base_flow_duthetadr * interpolated_RS *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    k * base_flow_duthetadz * interpolated_ZS *
                                    testf_ * Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_ur * dtestfdr * interpolated_RC *
                                    Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    base_flow_ur * dtestfdz * interpolated_ZC *
                                    Jbar * w / r;
            residuals[local_eqn] -= visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    interpolated_UC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * (1.0 + this->Gamma[0]) * k *
                                    interpolated_RC * base_flow_ur * testf_ *
                                    Jbar * w / (r * r);
            residuals[local_eqn] += k * visc_ratio * this->Gamma[0] *
                                    interpolated_WC * dtestfdz * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    base_flow_duzdr * interpolated_RC *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] -= k * visc_ratio * this->Gamma[0] *
                                    base_flow_duzdz * interpolated_ZC *
                                    dtestfdz * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * base_flow_duthetadz * dtestfdZS * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * r * interpolated_dVdZS * dtestfdz * Jbar * w;
            residuals[local_eqn] +=
              visc_ratio * r * base_flow_duthetadz * JhatS * dtestfdz * w;
            residuals[local_eqn] -= visc_ratio * interpolated_RS *
                                    base_flow_duthetadz * dtestfdz * Jbar * w;
            residuals[local_eqn] += visc_ratio * this->Gamma[0] *
                                    interpolated_dVdRS * testf_ * Jbar * w;
            residuals[local_eqn] -=
              visc_ratio * k * interpolated_UC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * base_flow_durdr *
                                    interpolated_RC * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * k * base_flow_durdz *
                                    interpolated_ZC * testf_ * Jbar * w / r;
            residuals[local_eqn] -=
              visc_ratio * interpolated_VS * testf_ * Jbar * w / r;
            residuals[local_eqn] += visc_ratio * base_flow_utheta *
                                    interpolated_RS * testf_ * Jbar * w /
                                    (r * r);
            residuals[local_eqn] -=
              scaled_re_st * r * base_flow_duthetadt * testf_ * JhatS * w;
            residuals[local_eqn] -=
              scaled_re * base_flow_utheta * base_flow_ur * testf_ * JhatS * w;
            residuals[local_eqn] += r * body_force[2] * testf_ * JhatS * w;
            residuals[local_eqn] +=
              scaled_re_inv_fr * r * G[2] * testf_ * JhatS * w;
            residuals[local_eqn] += k * base_flow_p * testf_ * JhatC * w;
            residuals[local_eqn] -= k * visc_ratio * (1.0 + this->Gamma[0]) *
                                    base_flow_ur * testf_ * JhatC * w / r;
            residuals[local_eqn] -=
              visc_ratio * base_flow_utheta * testf_ * JhatS * w / r;


            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions again
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[0]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * psif[l2] * dtestfdr *
                    Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * psif[l2] *
                    testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * k * psif[l2] * testf_ * Jbar * w / r;
                }

                // Radial velocity component (sine part) U_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duthetadr * testf_ *
                    Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * base_flow_utheta * psif[l2] * testf_ * Jbar * w;
                }

                // Axial velocity component (cosine part) W_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * visc_ratio * this->Gamma[0] * psif[l2] * dtestfdz *
                    Jbar * w;
                }

                // Axial velocity component (sine part) W_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[3]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * psif[l2] * base_flow_duthetadz * testf_ *
                    Jbar * w;
                }

                // Azimuthal velocity component (cosine part) V_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[4]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * scaled_re * base_flow_utheta * psif[l2] * testf_ *
                    Jbar * w;
                }

                // Azimuthal velocity component (sine part) V_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[5]);
                if (local_unknown >= 0)
                {
                  if (flag == 2)
                  {
                    // Add the mass matrix
                    mass_matrix(local_eqn, local_unknown) +=
                      scaled_re_st * r * psif[l2] * testf_ * Jbar * w;
                  }

                  // Add contributions to the Jacobian matrix
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_B[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_uz * group_A[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[1] * group_A[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_B[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * psif[l2] * dtestfdr * Jbar *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * k * psif[l2] *
                    testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_A[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * group_B[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * psif[l2] * testf_ * Jbar * w / r;
                }

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * base_flow_duthetadr *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * base_flow_durdr *
                    psif[l2] * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * base_flow_p * dtestfdr * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * base_flow_ur *
                    dtestfdr * psif[l2] * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * psif[l2] *
                    base_flow_ur * testf_ * Jbar * w / (r * r);
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * base_flow_duzdr *
                    psif[l2] * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * base_flow_durdr * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    k * base_flow_p * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_B[l2] * w / r;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * psif[l2] * base_flow_duthetadt * testf_ *
                    Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duthetadr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_ur * base_flow_duthetadr *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[0] *
                    base_flow_duthetadr * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re * r * base_flow_uz * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * mesh_velocity[1] * group_E[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * psif[l2] * base_flow_uz * base_flow_duthetadz *
                    testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * psif[l2] * mesh_velocity[1] *
                    base_flow_duthetadz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * body_force[2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * psif[l2] * G[2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadr * group_B[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_duthetadr *
                                                        dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * k *
                    base_flow_duthetadr * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadz * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * group_E[l2] * dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadz * group_B[l2] *
                    dtestfdz * w;
                  jacobian(local_eqn, local_unknown) -= visc_ratio * psif[l2] *
                                                        base_flow_duthetadz *
                                                        dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * base_flow_utheta * psif[l2] * testf_ * Jbar *
                    w / (r * r);
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duthetadt * testf_ *
                    group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * base_flow_utheta * base_flow_ur * testf_ *
                    group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[2] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[2] * testf_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * base_flow_utheta * testf_ * group_B[l2] * w /
                    r;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * scaled_re * base_flow_utheta * base_flow_duthetadz *
                    psif[l2] * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * base_flow_durdz *
                    psif[l2] * dtestfdr * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    k * base_flow_p * dtestfdz * psif[l2] * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * base_flow_ur *
                    dtestfdz * psif[l2] * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * this->Gamma[0] * base_flow_duzdz *
                    psif[l2] * dtestfdz * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * k * base_flow_durdz * psif[l2] * testf_ *
                    Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    k * base_flow_p * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    k * visc_ratio * (1.0 + this->Gamma[0]) * base_flow_ur *
                    testf_ * group_A[l2] * w / r;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * r * base_flow_ur * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * mesh_velocity[0] * group_E[l2] * testf_ *
                    w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_st * r * psif[l2] *
                    this->node_pt(l2)->position_time_stepper_pt()->weight(1,
                                                                          0) *
                    base_flow_duthetadz * testf_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * base_flow_duthetadr * group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * r * group_E[l2] * dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadr * group_A[l2] *
                    dtestfdr * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * base_flow_utheta *
                    group_F(l, l2) * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * (1.0 + this->Gamma[0]) * k * k *
                    base_flow_duthetadz * psif[l2] * testf_ * Jbar * w / r;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * r * base_flow_duthetadz * group_A[l2] *
                    dtestfdz * w;
                  jacobian(local_eqn, local_unknown) +=
                    visc_ratio * this->Gamma[0] * group_E[l2] * testf_ * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re_st * r * base_flow_duthetadt * testf_ *
                    group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    scaled_re * base_flow_utheta * base_flow_ur * testf_ *
                    group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    r * body_force[2] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) +=
                    scaled_re_inv_fr * r * G[2] * testf_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    visc_ratio * base_flow_utheta * testf_ * group_A[l2] * w /
                    r;
                }

              } // End of loop over velocity shape functions

              // Now loop over pressure shape functions
              // (This is the contribution from pressure gradient)
              for (unsigned l2 = 0; l2 < n_pres; l2++)
              {
                // Cosine part P_k^C
                local_unknown = this->p_local_eqn(l2, 0);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * psip[l2] * testf_ * Jbar * w;
                }

                // Sine part P_k^S has no contribution

              } // End of loop over pressure shape functions

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

        } // End of loop over fluid test functions


        // ====================
        // CONTINUITY EQUATIONS
        // ====================

        // Loop over the pressure test functions
        for (unsigned l = 0; l < n_pres; l++)
        {
          // Cache the test function
          const double testp_ = testp[l];

          // --------------------------------------
          // FIRST CONTINUITY EQUATION: COSINE PART
          // --------------------------------------

          // Get local equation number of first pressure value at this node
          local_eqn = this->p_local_eqn(l, 0);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += r * interpolated_dUdRC * testp_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RC * base_flow_durdr * testp_ * Jbar * w;
            residuals[local_eqn] += interpolated_UC * testp_ * Jbar * w;
            residuals[local_eqn] += k * interpolated_VS * testp_ * Jbar * w;
            residuals[local_eqn] -=
              k * base_flow_duthetadr * interpolated_RS * testp_ * Jbar * w;
            residuals[local_eqn] -=
              k * base_flow_duthetadz * interpolated_ZS * testp_ * Jbar * w;
            residuals[local_eqn] += r * interpolated_dWdZC * testp_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RC * base_flow_duzdz * testp_ * Jbar * w;
            residuals[local_eqn] -=
              interpolated_RC * source * testp_ * Jbar * w;
            residuals[local_eqn] += base_flow_ur * testp_ * JhatC * w;
            residuals[local_eqn] -= source * r * testp_ * JhatC * w;


            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[0]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * group_B[l2] * testp_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * testp_ * Jbar * w;
                }

                // Radial velocity component (sine part) U_k^S
                // has no contribution

                // Axial velocity component (cosine part) W_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[2]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * group_A[l2] * testp_ * w;
                }

                // Axial velocity component (sine part) W_k^S
                // has no contribution

                // Azimuthal velocity component (cosine part) V_k^C
                // has no contribution

                // Azimuthal velocity component (sine part) V_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[5]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * psif[l2] * testp_ * Jbar * w;
                }

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_durdr * testp_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    r * group_D[l2] * testp_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_duzdz * testp_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    psif[l2] * source * testp_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_ur * testp_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    source * r * testp_ * group_B[l2] * w;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * base_flow_duthetadr * psif[l2] * testp_ * Jbar * w;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * group_C[l2] * testp_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_ur * testp_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    source * r * testp_ * group_A[l2] * w;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * base_flow_duthetadz * psif[l2] * testp_ * Jbar * w;
                }

              } // End of loop over velocity shape functions

              // Real and imaginary pressure components, P_k^C and P_k^S,
              // have no contribution

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

          // -------------------------------------
          // SECOND CONTINUITY EQUATION: SINE PART
          // -------------------------------------

          // Get local equation number of second pressure value at this node
          local_eqn = this->p_local_eqn(l, 1);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += r * interpolated_dUdRS * testp_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RS * base_flow_durdr * testp_ * Jbar * w;
            residuals[local_eqn] += interpolated_US * testp_ * Jbar * w;
            residuals[local_eqn] -= k * interpolated_VC * testp_ * Jbar * w;
            residuals[local_eqn] +=
              k * base_flow_duthetadr * interpolated_RC * testp_ * Jbar * w;
            residuals[local_eqn] +=
              k * base_flow_duthetadz * interpolated_ZC * testp_ * Jbar * w;
            residuals[local_eqn] += r * interpolated_dWdZS * testp_ * Jbar * w;
            residuals[local_eqn] +=
              interpolated_RS * base_flow_duzdz * testp_ * Jbar * w;
            residuals[local_eqn] -=
              interpolated_RS * source * testp_ * Jbar * w;
            residuals[local_eqn] += base_flow_ur * testp_ * JhatS * w;
            residuals[local_eqn] -= source * r * testp_ * JhatS * w;


            // Calculate the Jacobian
            // ----------------------

            if (flag)
            {
              // Loop over the velocity shape functions
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Radial velocity component (cosine part) U_k^C
                // has no contribution

                // Radial velocity component (sine part) U_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[1]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * group_B[l2] * testp_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * testp_ * Jbar * w;
                }

                // Axial velocity component (cosine part) W_k^C
                // has no contribution

                // Axial velocity component (sine part) W_k^S
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[3]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * group_A[l2] * testp_ * w;
                }

                // Azimuthal velocity component (cosine part) V_k^C
                local_unknown = this->nodal_local_eqn(l2, u_nodal_index[4]);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) -=
                    k * psif[l2] * testp_ * Jbar * w;
                }

                // Azimuthal velocity component (sine part) V_k^S
                // has no contribution

                // Perturbation to radial nodal coord (cosine part) R_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 0));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * base_flow_duthetadr * psif[l2] * testp_ * Jbar * w;
                }

                // Perturbation to radial nodal coord (sine part) R_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 1));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_durdr * testp_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    r * group_D[l2] * testp_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * base_flow_duzdz * testp_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) -=
                    psif[l2] * source * testp_ * Jbar * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_ur * testp_ * group_B[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    source * r * testp_ * group_B[l2] * w;
                }

                // Perturbation to axial nodal coord (cosine part) Z_k^C
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 2));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    k * base_flow_duthetadz * psif[l2] * testp_ * Jbar * w;
                }

                // Perturbation to axial nodal coord (sine part) Z_k^S
                local_unknown = this->nodal_local_eqn(
                  l2, this->xhat_index_lin_axi_nst(l2, 3));
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    r * group_C[l2] * testp_ * w;
                  jacobian(local_eqn, local_unknown) +=
                    base_flow_ur * testp_ * group_A[l2] * w;
                  jacobian(local_eqn, local_unknown) -=
                    source * r * testp_ * group_A[l2] * w;
                }

              } // End of loop over velocity shape functions

              // Real and imaginary pressure components, P_k^C and P_k^S,
              // have no contribution

              // Geometric contribution to jacobian
              // TODO

            } // End of Jacobian calculation

          } // End of if not boundary condition statement

        } // End of loop over pressure test functions

        // TOTAL VELOCITY EQUATIONS
        //-------------------

        // Loop over the velocity test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < this->n_u_lin_axi_nst(); i++)
          {
            // Additional velocity data
            // ------------------------

            // Find its local equation number
            local_eqn = this->nodal_local_eqn(l, u_index_lin_axi_nst_fe(l, i));

            // If it is not pinned
            if (local_eqn >= 0)
            {
              Vector<double> pos_n(2, 0.0);
              for (unsigned k = 0; k < 2; k++)
              {
                pos_n[k] = this->nodal_position(l, k);
              }

              residuals[local_eqn] +=
                (this->nodal_value(l, this->u_index_lin_axi_nst(i)) -
                 (this->nodal_value(l, u_index_lin_axi_nst_fe(l, i)) +
                  u_bar(pos_n, i)));
            }
          } // End of loop over velocity components
        } // End of loop over test functions

        // TOTAL PRESSURE EQUATION
        //-------------------

        // Loop over the Nodes
        for (unsigned l = 0; l < this->npres_lin_axi_nst(); l++)
        {
          for (unsigned j = 0; j < 2; j++)
          {
            // Get the local equation number
            local_eqn = this->nodal_local_eqn(l, p_index_lin_axi_nst_fe(l, j));

            // If not a boundary conditions
            if (local_eqn >= 0)
            {
              // If not subject to Dirichlet BC
              Vector<double> pos_n(2, 0.0);
              for (unsigned k = 0; k < 2; k++)
              {
                pos_n[k] = this->nodal_position(l, k);
              }

              residuals[local_eqn] +=
                (this->nodal_value(l, this->p_index_lin_axi_nst(j)) -
                 (this->nodal_value(l, p_index_lin_axi_nst_fe(l, j)) + p_bar(pos_n)));
            }
          }
        } // End of loop over l
      } // End of loop over the integration points
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
        interpolated_u +=
          this->nodal_value(l, u_index_lin_axi_nst_fe(l, i)) * psi[l];
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
        interpolated_p +=
          this->nodal_value(l, p_index_lin_axi_nst_fe(l, i)) * psi[l];
      }

      return (interpolated_p);
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
        Vector<double> x(2, 0.0);
        this->interpolated_x(s, x);

        // Output global coordinates to file
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << x[i] << " ";
        }

        // Output perturbations to nodal positions to file
        for (unsigned i = 0; i < 4; i++)
        {
          outfile << this->interpolated_nodal_position_perturbation_lin_axi_nst(
                       s, i)
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
          outfile << u_bar(x, i) << " ";
        }

        // Output pressure to file
        double p = p_bar(x);
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << p << " ";
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
