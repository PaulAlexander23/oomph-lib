using namespace std;
using namespace oomph;

namespace oomph
{
  //=======================================================================
  /// 1D Free surface elements for Hele Shaw problems with pseudo-elastic
  /// mesh motion
  //======================================================================
  template<class ELEMENT>
  class HeleShawInterfaceElement : public virtual FaceGeometry<ELEMENT>,
                                   public virtual FaceElement
  {
  public:
    typedef void (*DropPressureFctPt)(const Vector<double>& x,
                                      double& p_bubble);

    typedef void (*DropPressureGradientFctPt)(const Vector<double>& x,
                                              Vector<double>& grad_pb);

    /// \short Function pointer to function which provides h(x,t) and dh/dt, for
    /// vector x.
    typedef void (*UpperWallFctPt)(const Vector<double>& x,
                                   double& h,
                                   double& dhdt);

    typedef void (*UpperWallFluxFctPt)(const Vector<double>& x,
                                       double& h,
                                       double& dhdt,
                                       Vector<double>& dhdx,
                                       Vector<double>& d_dhdt_dx);

    /// \short Function pointer to function which provides wall speed as a
    /// function of x. This allows us to solve for bubble motion in a moving
    /// frame. A constant wall speed does not affect the mass conservation
    /// equations, but does feature in the kinematic equation for interface
    /// motion.
    typedef void (*WallSpeedFctPt)(const Vector<double>& x,
                                   Vector<double>& U_wall);

    WallSpeedFctPt& wall_speed_fct_pt()
    {
      return Wall_speed_fct_pt;
    }
    WallSpeedFctPt wall_speed_fct_pt() const
    {
      return Wall_speed_fct_pt;
    }

    DropPressureFctPt& drop_pressure_fct_pt()
    {
      return Drop_pressure_fct_pt;
    }
    DropPressureFctPt drop_pressure_fct_pt() const
    {
      return Drop_pressure_fct_pt;
    }
    DropPressureGradientFctPt& drop_pressure_gradient_fct_pt()
    {
      return Drop_pressure_gradient_fct_pt;
    }

    //    unsigned which_bubble;
  private:
    WallSpeedFctPt Wall_speed_fct_pt;
    DropPressureFctPt Drop_pressure_fct_pt;
    DropPressureGradientFctPt Drop_pressure_gradient_fct_pt;

    /// \short Pointer to the alpha ratio:
    double* Alpha_pt;

    double* Drop_viscosity_pt;

    /// Pointer to Q^-1
    double* Q_inv_pt;

    /// Pointer to the Strouhal number
    double* St_pt;

    double* Measure_pt;


    /// Default value for physical constants
    static double Default_Physical_Constant_Value;

    /// \short ID of additional unknowns introduced by this face element
    /// (smoothed components of derivative of tangent vector, and
    /// Lagrange multiplier)
    unsigned Id;

    /// \short Index of the nodal value at which the pressure is stored
    unsigned P_index_bulk;
    unsigned P_index_drop;

    /// \short Equation number of the equation for the Lagrange multiplier
    int lagrange_local_eqn(const unsigned& j)
    {
      // Get the index of the nodal value associated with Lagrange multiplier
      // NOTE: It's the first new one
      unsigned lagr_index =
        dynamic_cast<BoundaryNodeBase*>(node_pt(j))
          ->index_of_first_value_assigned_by_face_element(Id) +
        0;

      // Return nodal value
      return this->nodal_local_eqn(j, lagr_index);
    }

    /// Return the Lagrange multiplier at local node j
    double& lagrange(const unsigned& j)
    {
      // Get the index of the nodal value associated with Lagrange multiplier
      // NOTE: It's the first new one
      unsigned lagr_index =
        dynamic_cast<BoundaryNodeBase*>(node_pt(j))
          ->index_of_first_value_assigned_by_face_element(Id) +
        0;

      // Return (ref) to value
      return *node_pt(j)->value_pt(lagr_index);
    }

    /// \short Equation number of equation that does the projection
    /// for the derivative of the i-th component of the tangent
    /// vector at node j
    int projected_tangent_deriv_local_eqn(const unsigned& j, const unsigned& i)
    {
      // Get the index of the nodal value associated with projected tang vector
      unsigned tang_index =
        dynamic_cast<BoundaryNodeBase*>(node_pt(j))
          ->index_of_first_value_assigned_by_face_element(Id) +
        1 + i;

      // Return nodal value
      return this->nodal_local_eqn(j, tang_index);
    }


    /// \short Return i-th component of projected deriv of tangent vector
    /// at local node j
    double& projected_tangent_deriv(const unsigned& j, const unsigned& i)
    {
      // Get the index of the nodal value associated with tang deriv
      unsigned veloc_index =
        dynamic_cast<BoundaryNodeBase*>(node_pt(j))
          ->index_of_first_value_assigned_by_face_element(Id) +
        1 + i;

      // Return ref to value
      return *node_pt(j)->value_pt(veloc_index);
    }

    /// \short Helper function to calculate the residuals and
    /// (if flag==1) the Jacobian of the equations.
    void fill_in_generic_residual_contribution_hele_shaw_interface(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag);


  public:
    /// \short Constructor, pass a pointer to the bulk element and the face
    /// index of the bulk element to which the element is to be attached to. The
    /// optional identifier can be used to distinguish the additional nodal
    /// values created by this element (in order: 0: Lagrange multiplier; 1:
    /// x-component of smoothed derivative of tangent vector; 2: y-component of
    /// smoothed derivative of tangent vector)
    /// from those created by other FaceElements.
    HeleShawInterfaceElement(FiniteElement* const& element_pt,
                             const int& face_index,
                             const unsigned& id = 0)
      :

        FaceGeometry<ELEMENT>(),
        Id(id)
    {
      // Attach the geometrical information to the element
      // This function also assigned nbulk_value from required_nvalue of the
      // bulk element
      element_pt->build_face_element(face_index, this);

      Wall_speed_fct_pt = 0;
      Drop_pressure_fct_pt = 0;

      P_index_bulk = 0;
      P_index_drop = 1;

      // Read out the number of nodes on the face
      unsigned n_node_face = this->nnode();

      // Set the additional data values in the face
      // There are three additional values at each node -- the
      // Lagrange multiplier and the two components of the smoothed
      // derivative of the tangent vector
      /// Now we add another lagrange multiplier for tangential stress.
      Vector<unsigned> additional_data_values(n_node_face);
      for (unsigned i = 0; i < n_node_face; i++)
      {
        additional_data_values[i] = 3;
      }

      // Now add storage and set the map containing
      // the position of the first entry of this face element's
      // additional values.
      add_additional_values(additional_data_values, id);

      //        which_bubble=0;
    }

    /// Access function: Pointer to source function
    UpperWallFctPt& upper_wall_fct_pt()
    {
      return Upper_wall_fct_pt;
    }

    /// Access function: Pointer to source function. Const version
    UpperWallFctPt upper_wall_fct_pt() const
    {
      return Upper_wall_fct_pt;
    }

    UpperWallFluxFctPt& upper_wall_flux_fct_pt()
    {
      return Upper_wall_flux_fct_pt;
    }

    UpperWallFluxFctPt upper_wall_flux_fct_pt() const
    {
      return Upper_wall_flux_fct_pt;
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // std::cout << "Interface"<< std::endl;
      //        /// Fill in the Jacobian contribution
      //        fill_in_generic_residual_contribution_hele_shaw_interface(residuals,jacobian,mass_matrix,1);
      fill_in_contribution_to_jacobian(residuals, jacobian);

      /// Fill in the mass matrix contribution
      fill_in_generic_residual_contribution_hele_shaw_interface(
        residuals, jacobian, mass_matrix, 2);
      //        /// Alice: For my sanity later on. External dependencies must
      //        still be added in driver code.
      //        /// May be useful for finding steady solutions where frame speed
      //        or capillary number
      //        /// are free parameters. The residuals depend linearly on these
      //        parameters, so fd
      //        /// is not completely unreasonable.
      //        /// NB. Similar hack in Poisson flux elements?
      //        this->fill_in_jacobian_from_external_by_fd(jacobian);
    }


    /// Calculate the residuals by calling the generic residual contribution.
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Add the residual contributions
      fill_in_generic_residual_contribution_hele_shaw_interface(
        residuals,
        GeneralisedElement::Dummy_matrix,
        GeneralisedElement::Dummy_matrix,
        0);
    }


    /// \short The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    /// \short Virtual function that specifies the non-dimensional
    /// surface tension as a function of local position within the element.
    /// The default behaviour is a constant surface tension of value 1.0
    /// This function can be overloaded in more
    /// specialised elements to incorporate variations in surface tension.
    virtual double sigma(const Vector<double>& s_local)
    {
      return 1.0;
    }

    /// The value of the inverse Capillary number
    const double& q_inv() const
    {
      return *Q_inv_pt;
    }

    /// Pointer to the inverse Capillary number
    double*& q_inv_pt()
    {
      return Q_inv_pt;
    }

    /// The value of the Strouhal number
    const double& st() const
    {
      return *St_pt;
    }

    /// The pointer to the Strouhal number
    double*& st_pt()
    {
      return St_pt;
    }

    double*& measure_pt()
    {
      return Measure_pt;
    }

    double alpha() const
    {
      if (Alpha_pt == 0)
      {
        return 1.0;
      }
      return *Alpha_pt;
    }

    double*& alpha_pt()
    {
      return Alpha_pt;
    }
    double* alpha_pt() const
    {
      return Alpha_pt;
    }

    /// Return the value of the external pressure
    double get_p_bubble(Vector<double>& x) const
    {
      if (Drop_pressure_fct_pt == 0)
      {
        return 0.0;
      }
      else
      {
        double pressure;
        (*Drop_pressure_fct_pt)(x, pressure);
        return pressure;
      }
    }

    double drop_viscosity() const
    {
      if (Drop_viscosity_pt == 0)
      {
        return 1.0;
      }
      return *Drop_viscosity_pt;
    }


    double*& drop_viscosity_pt()
    {
      return Drop_viscosity_pt;
    }
    double* drop_viscosity_pt() const
    {
      return Drop_viscosity_pt;
    }

    void get_wall_velocity(Vector<double>& x, Vector<double>& U)
    {
      if (Wall_speed_fct_pt == 0)
      {
        U[0] = 0.0;
        U[1] = 0.0;
      }
      else
      {
        (*Wall_speed_fct_pt)(x, U);
      }
    }

    /// Overload the output functions
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// Output the element
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Local coordinate
      Vector<double> s(1);

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Get the value of the Capillary number
      double Q_inv = q_inv();

      outfile << "ZONE\n";

      // Loop over plot points
      for (unsigned i_plot = 0; i_plot < n_plot; i_plot++)
      {
        // Get coordinate
        get_s_plot(i_plot, n_plot, s);


        // Call the derivatives of the shape function at the knot point
        this->dshape_local(s, psif, dpsifds);

        // Compute what we need...
        Vector<double> interpolated_x(2, 0.0);
        Vector<double> interpolated_dx_dt(2, 0.0);
        Vector<double> tau(2, 0.0);

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            // Smoothed tangent vector
            tau[i] += projected_tangent_deriv(l, i) * psif(l);

            // Spatial bits
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_dx_dt[i] += this->dnodal_position_dt(l, i) * psif(l);
          }
        }


        // Now calculate the unit normal vector
        Vector<double> interpolated_n(2);
        outer_unit_normal(s, interpolated_n);

        // Also get the (possibly variable) surface tension
        double sigma_local = this->sigma(s);

        // Assemble the surface tension and normal speed terms
        double sigma_kappa = 0.0;
        for (unsigned k = 0; k < 2; k++)
        {
          sigma_kappa += sigma_local * Q_inv * tau[k] * interpolated_n[k];
        }

        outfile << interpolated_x[0] << " " << interpolated_x[1] << " "
                << tau[0] << " " << tau[1] << " " << interpolated_dx_dt[0]
                << " " << interpolated_dx_dt[1] << " " << sigma_kappa << " "
                << std::endl;
      }
    }

    /// Overload the C-style output function
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// C-style Output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      // hierher fill this in when we've finalised what lives in
      // output
    }

    /// Pointer to function that specifies the gap width and wall velocity
    UpperWallFctPt Upper_wall_fct_pt;

    UpperWallFluxFctPt Upper_wall_flux_fct_pt;
  };


  //=======================================================================
  /// Calculate the residual contribution (kinematic and dynamic BC and
  /// Lagrange multiplier contribution from pseudo-elastic node updates
  //=======================================================================
  template<class ELEMENT>
  void HeleShawInterfaceElement<ELEMENT>::
    fill_in_generic_residual_contribution_hele_shaw_interface(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // std::cout << "New element " << std::endl;
    if (flag == 1)
    {
      std::cout << "Flag is 1!" << std::endl;
    }

    int measure_eqn = external_local_eqn(0, 0);

    if (0)
    {
      std::cout << node_pt(0)->value(0) << " " << node_pt(0)->value(1) << " "
                << node_pt(0)->value(2) << " " << node_pt(0)->value(3) << " "
                << node_pt(0)->value(4) << " " << std::endl;
    }

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, 1);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Get the value of the Capillary number
    double Q_inv = q_inv();

    // Get the value of the Strouhal numer
    double St = st();

    // Integers to store the local equation numbers
    int local_eqn = 0;
    int local_eqn_drop;
    int local_eqn_bulk;


    // Storage for the local cooridinate
    Vector<double> s(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the local coordinate at the integration point
      s[0] = integral_pt()->knot(ipt, 0);

      // Get the integral weight
      double W = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape function at the knot point
      this->dshape_local_at_knot(ipt, psif, dpsifds);

      // Compute what we need...
      double interpolated_p_bulk = 0.0;
      double interpolated_p_drop = 0.0;
      double interpolated_lagrange = 0.0;
      Vector<double> interpolated_tangent(2, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> interpolated_dx_dt(2, 0.0);
      Vector<double> tau(2, 0.0);

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_lagrange += lagrange(l) * psif[l];
        interpolated_p_bulk += node_pt(l)->value(P_index_bulk) * psif[l];
        interpolated_p_drop += node_pt(l)->value(P_index_drop) * psif[l];

        // Loop over directional components
        for (unsigned i = 0; i < 2; i++)
        {
          // Smoothed tangent vector
          tau[i] += projected_tangent_deriv(l, i) * psif(l);

          // Spatial bits
          interpolated_x[i] += this->nodal_position(l, i) * psif(l);
          interpolated_dx_dt[i] += this->dnodal_position_dt(l, i) * psif(l);
          interpolated_tangent[i] += this->nodal_position(l, i) * dpsifds(l, 0);
        }
      }

      // Get the value of the bubble pressure. This need not be constant, and in
      // fact currently
      // contains the transverse curvature term.
      double p_bubble = get_p_bubble(interpolated_x);

      Vector<double> U_wall(2, 0.0);
      get_wall_velocity(interpolated_x, U_wall);

      // Calculate the length of the tangent Vector
      double tlength = interpolated_tangent[0] * interpolated_tangent[0] +
                       interpolated_tangent[1] * interpolated_tangent[1];

      // Set the Jacobian of the line element
      double J = sqrt(tlength);

      // Normalise the tangent Vector
      interpolated_tangent[0] /= J;
      interpolated_tangent[1] /= J;
      double test_sign = interpolated_tangent[0] * interpolated_x[1] -
                         interpolated_tangent[1] * interpolated_x[0];
      //	std::cout <<
      //interpolated_tangent[0]*interpolated_x[1]-interpolated_tangent[1]*interpolated_x[0]
      //<< std::endl;
      if (test_sign < 0)
      {
        //	interpolated_tangent[0]*=-1;
        //	interpolated_tangent[1]*=-1;
      }

      // Now calculate the unit normal vector
      Vector<double> interpolated_n(2);
      outer_unit_normal(ipt, interpolated_n);

      //	double n_dot_x =
      //interpolated_n[0]*interpolated_x[0]+interpolated_n[1]*interpolated_x[1];
      //	std::cout << n_dot_x << std::endl;

      double kappa = 0.0;
      double normal_speed_time = 0.0;
      double normal_speed_wall = 0.0;
      for (unsigned k = 0; k < 2; k++)
      {
        //     sigma_kappa+=sigma_local*Ca_inv*tau[k]*interpolated_n[k];
        kappa += tau[k] * interpolated_n[k];
        normal_speed_time += interpolated_dx_dt[k] * interpolated_n[k];
        normal_speed_wall += U_wall[k] * interpolated_n[k];
      }

      double h = 0.0;
      double dhdt = 0.0;
      upper_wall_fct_pt()(interpolated_x, h, dhdt);

      // Non-dim gap-width
      double local_alpha = alpha();
      double local_drop_viscosity = drop_viscosity();

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Eqns to determine the smoothed derivatives of the tangent vector
        // ----------------------------------------------------------------
        for (unsigned i = 0; i < 2; i++)
        {
          local_eqn = projected_tangent_deriv_local_eqn(l, i);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              (tau[i] * psif(l) * J + interpolated_tangent[i] * dpsifds(l, 0)) *
              W;
          }
        }

        // Eqn for Lagrange multiplier (dynamic free surface condition)
        //-------------------------------------------------------------
        local_eqn = lagrange_local_eqn(l);
        if (local_eqn >= 0)
        {
          double Q = 1 / Q_inv;
          double alpha = local_alpha;
          residuals[local_eqn] +=
            (interpolated_p_bulk - interpolated_p_drop - p_bubble +
             kappa / (3.0 * alpha * alpha * Q)) *
            psif(l) * W * J;
        }

        // Contribution to arbitrary solution measure
        //---------------------------------------------
        if (measure_eqn >= 0)
        {
          std::cout << "Measure eqn" << std::endl;
          residuals[measure_eqn] += 0;
        }


        // Contribution to bulk equation (kinematic bc)
        //---------------------------------------------
        local_eqn_bulk = nodal_local_eqn(l, P_index_bulk);
        local_eqn_drop = nodal_local_eqn(l, P_index_drop);
        double mu_drop = local_drop_viscosity;
        if (local_eqn_bulk >= 0)
        {
          residuals[local_eqn_bulk] +=
            h * St * normal_speed_time * psif(l) * W * J;
          residuals[local_eqn_bulk] += h * normal_speed_wall * psif(l) * W * J;
          residuals[local_eqn_drop] -=
            mu_drop * h * St * normal_speed_time * psif(l) * W * J;
          residuals[local_eqn_drop] -=
            mu_drop * h * normal_speed_wall * psif(l) * W * J;


          /// These residuals depend only on solid position, and
          /// possibly external data. (but not nodal data)

          /// True mass matrix.
          if (flag == 2)
          {
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              /// Do x-components.
              unsigned local_unknown = this->position_local_eqn(l2, 0, 0);
              if (local_unknown >= 0)
              {
                mass_matrix(local_eqn_bulk, local_unknown) +=
                  h * St * psif(l) * W * J * interpolated_n[0] * psif(l2);
                mass_matrix(local_eqn_drop, local_unknown) -=
                  mu_drop * St * psif(l) * W * J * interpolated_n[0] * psif(l2);
              }

              /// Do y-components
              local_unknown = this->position_local_eqn(l2, 0, 1);
              if (local_unknown >= 0)
              {
                mass_matrix(local_eqn_bulk, local_unknown) +=
                  h * St * psif(l) * W * J * interpolated_n[1] * psif(l2);
                mass_matrix(local_eqn_drop, local_unknown) -=
                  mu_drop * h * St * psif(l) * W * J * interpolated_n[1] *
                  psif(l2);
              }
            }
          }
        }

        // Lagrange multiplier contributions to pseudo_solid equations
        //------------------------------------------------------------
        for (unsigned i = 0; i < 2; i++)
        {
          local_eqn = this->position_local_eqn(l, 0, i);

          if (local_eqn >= 0)
          {
            residuals[local_eqn] +=
              -interpolated_lagrange * interpolated_n[i] * psif[l] * W * J;
          }
        }

      } // End of loop over shape functions
    } // End of loop over integration points

    if (0)
    {
      unsigned l = 1;
      int local_eqn = 0;
      for (unsigned i = 0; i < 2; i++)
      {
        int local_eqn = projected_tangent_deriv_local_eqn(l, i);
        std::cout << residuals[local_eqn] << " ";
      }
      /// Dynamic boundary condition
      local_eqn = lagrange_local_eqn(l);
      std::cout << residuals[local_eqn] << " ";

      local_eqn = nodal_local_eqn(l, P_index_bulk);
      std::cout << residuals[local_eqn] << " ";
      local_eqn = nodal_local_eqn(l, P_index_drop);
      std::cout << residuals[local_eqn] << " ";
      std::cout << std::endl;
    }
  }


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


} // namespace oomph
