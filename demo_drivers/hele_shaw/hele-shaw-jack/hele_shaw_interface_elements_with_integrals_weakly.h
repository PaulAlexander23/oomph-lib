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
    /// \short Function pointer to function which provides bubble pressure as a
    /// function of vector x. This should usually be a constant - should default
    /// to zero. Any perturbations to the system can be supplied through a
    /// spatially varying bubble pressure function.
    typedef void (*BubblePressureFctPt)(const Vector<double>& x,
                                        double& p_bubble);

    typedef void (*BubblePressureGradientFctPt)(const Vector<double>& x,
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


    BubblePressureFctPt& bubble_pressure_fct_pt()
    {
      return Bubble_pressure_fct_pt;
    }
    BubblePressureFctPt bubble_pressure_fct_pt() const
    {
      return Bubble_pressure_fct_pt;
    }
    BubblePressureGradientFctPt& bubble_pressure_gradient_fct_pt()
    {
      return Bubble_pressure_gradient_fct_pt;
    }

    WallSpeedFctPt& wall_speed_fct_pt()
    {
      return Wall_speed_fct_pt;
    }
    WallSpeedFctPt wall_speed_fct_pt() const
    {
      return Wall_speed_fct_pt;
    }

    bool* construct_mass_matrix_left_multiplier_pt;

    bool* use_fd_for_dresidual_dnodal_coordinates_pt;

    bool* fill_solid_by_fd_pt;

  private:
    /// Pointer to function that specifies the gap width and wall velocity
    BubblePressureFctPt Bubble_pressure_fct_pt;
    BubblePressureGradientFctPt Bubble_pressure_gradient_fct_pt;

    WallSpeedFctPt Wall_speed_fct_pt;

    /// \short Pointer to the alpha ratio: h/L
    double* Alpha_pt;

    /// Pointer to Q^-1
    double* Q_inv_pt;

    /// Pointer to Q_zero
    int* Q_zero_pt;

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
    unsigned P_index_interface;

    void get_dresidual_dnodal_coordinates(
      RankThreeTensor<double>& dresidual_dnodal_coordinates);

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

    /// \short Equation number of the equation for the Lagrange multiplier
    int tangential_lagrange_local_eqn(const unsigned& j)
    {
      // Get the index of the nodal value associated with Lagrange multiplier
      // NOTE: It's the first new one
      unsigned tangential_lagr_index =
        dynamic_cast<BoundaryNodeBase*>(node_pt(j))
          ->index_of_first_value_assigned_by_face_element(Id) +
        3;


      // Return nodal value
      return this->nodal_local_eqn(j, tangential_lagr_index);
    }

    /// Return the Lagrange multiplier at local node j
    double& tangential_lagrange(const unsigned& j)
    {
      // Get the index of the nodal value associated with Lagrange multiplier
      // NOTE: It's the first new one
      unsigned tangential_lagr_index =
        dynamic_cast<BoundaryNodeBase*>(node_pt(j))
          ->index_of_first_value_assigned_by_face_element(Id) +
        3;

      //        double value = *node_pt(j)->value_pt(tangential_lagr_index);
      //        std::cout << value << std::endl;
      // Return (ref) to value
      return *node_pt(j)->value_pt(tangential_lagr_index);
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
      Bubble_pressure_fct_pt = 0;

      P_index_interface = 0;

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
        additional_data_values[i] = 4;
      }

      // Now add storage and set the map containing
      // the position of the first entry of this face element's
      // additional values.
      add_additional_values(additional_data_values, id);

      integral_measure_data_number = 0;
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


    unsigned integral_measure_data_number;


    /// Calculate the residuals by calling the generic residual contribution.
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Add the residual contributions
      fill_in_generic_residual_contribution_hele_shaw_interface(
        residuals,
        GeneralisedElement::Dummy_matrix,
        GeneralisedElement::Dummy_matrix,
        0); // JACK
    }

    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Add the residual contributions
      fill_in_generic_residual_contribution_hele_shaw_interface(
        residuals, jacobian, GeneralisedElement::Dummy_matrix, 1);
      this->fill_in_jacobian_from_external_by_fd(jacobian);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
      fill_in_generic_residual_contribution_hele_shaw_interface(
        residuals, jacobian, mass_matrix, 2);
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

    /// The value of the Q_zero
    const int& q_zero() const
    {
      return *Q_zero_pt;
    }

    /// Pointer to the inverse Capillary number
    double*& q_inv_pt()
    {
      return Q_inv_pt;
    }

    /// Pointer to Q_zero_pt
    int*& q_zero_pt()
    {
      return Q_zero_pt;
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


    /// \short Non-dimensional gap width (modulated by user-specifiable
    /// upper wall data)
    // hierher maybe get this from the bulk elements for all of these!
    double alpha() const
    {
      if (Alpha_pt == 0)
      {
        return 1.0;
      }
      return *Alpha_pt;
    }

    /// \short Pointer to non-dimensional gap width (modulated by
    /// user-specifiable upper wall data)
    double*& alpha_pt()
    {
      return Alpha_pt;
    }


    /// \short Pointer to non-dimensional gap width (modulated by
    /// user-specifiable upper wall data). Const version.
    double* alpha_pt() const
    {
      return Alpha_pt;
    }

    /// Return the value of the external pressure
    double get_p_bubble(Vector<double>& x) const
    {
      if (Bubble_pressure_fct_pt == 0)
      {
        return 0.0;
      }
      else
      {
        double pressure;
        (*Bubble_pressure_fct_pt)(x, pressure);
        return pressure;
      }
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
    int measure_eqn = external_local_eqn(0, 0);

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, 1);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Get the value of the Capillary number
    double Q_inv = q_inv();

    // Get the value of Q_zero
    int Q_zero = q_zero();

    // Get the value of the Strouhal numer
    double St = st();

    // Integers to store the local equation numbers
    int local_eqn = 0;
    int local_unknown;

    //    std::cout << "Measure is data " << integral_measure_data_number << "
    //    with " << external_data_pt(integral_measure_data_number)->nvalue() <<
    //    " values " << std::endl;

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
      double interpolated_p = 0.0;
      double interpolated_lagrange = 0.0;
      double interpolated_tangential_lagrange = 0.0;
      Vector<double> interpolated_tangent(2, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> interpolated_dx_dt(2, 0.0);
      //        double interpolated_Ys = 0.0;
      Vector<double> tau(2, 0.0);

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_lagrange += lagrange(l) * psif[l];
        interpolated_tangential_lagrange += tangential_lagrange(l) * psif[l];
        interpolated_p += node_pt(l)->value(P_index_interface) * psif[l];
        //            interpolated_Ys += node_pt(l)->value(4)*psif[l];

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
      // contains the transverse curvature term - JACK - this is where ther
      // transverse curvature term is contained, i.e.
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

      // Now calculate the unit normal vector
      Vector<double> interpolated_n(2);
      outer_unit_normal(ipt, interpolated_n);

      double kappa = 0.0;
      double normal_speed_time = 0.0;
      double normal_speed_wall = 0.0;
      double tangential_speed = 0.0;
      for (unsigned k = 0; k < 2; k++)
      {
        //     sigma_kappa+=sigma_local*Ca_inv*tau[k]*interpolated_n[k];
        kappa += tau[k] * interpolated_n[k];
        normal_speed_time += interpolated_dx_dt[k] * interpolated_n[k];
        normal_speed_wall += U_wall[k] * interpolated_n[k];
        tangential_speed += interpolated_dx_dt[k] * interpolated_tangent[k];
      }
      double direction = interpolated_n[0] * interpolated_tangent[1] -
                         interpolated_n[1] * interpolated_tangent[0];
      double sign = 0;
      if (direction > 0)
      {
        sign = 1;
      }
      if (direction < 0)
      {
        sign = -1;
      }

      double h = 0.0;
      double dhdt = 0.0;
      upper_wall_fct_pt()(interpolated_x, h, dhdt);

      // Non-dim gap-width
      double local_alpha = alpha();

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

            // Do the Jacobian calculation
            if (flag == 1)
            {
              // Loop over the nodes
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Derivatives w.r.t. solid positions will be handled by FDing.
                // Only 'tau' above is not solid position.
                local_unknown = projected_tangent_deriv_local_eqn(l2, i);
                if (local_unknown >= 0)
                {
                  //                            std::cout << "Interface 601" <<
                  //                            local_eqn << " " <<
                  //                            local_unknown << std::endl;
                  jacobian(local_eqn, local_unknown) +=
                    psif(l2) * psif(l) * W * J;
                }
              }
            } // End of Jacobian calculation
          }
        }

        // Eqn for Lagrange multiplier (dynamic free surface condition)
        //-------------------------------------------------------------
        local_eqn = lagrange_local_eqn(l);
        if (local_eqn >= 0)
        {
          double Q = 1 / Q_inv;
          double alpha = local_alpha;

          if (Q_zero == 1)
          {
            // Full equations - analysis in Q_inv
            residuals[local_eqn] +=
              (interpolated_p - p_bubble + kappa / (3.0 * alpha * alpha * Q)) *
              psif(l) * W * J;

            // Full equations - analysis in Q
            //    residuals[local_eqn] +=
            //(Q*(interpolated_p - p_bubble) +
            //kappa/(3.0*alpha*alpha))*psif(l)*W*J;

            //                         +
            //                         local_alpha*Q_inv*(local_alpha*kappa)/3.0)*psif(l)*W*J;
          }

          if (Q_zero == 0)
          {
            // Set Q_inv = 0
            residuals[local_eqn] +=
              (interpolated_p - p_bubble) * psif(l) * W * J;

            // Set Q = 0
            // residuals[local_eqn] += (kappa/(3.0*alpha*alpha))*psif(l)*W*J;
          }

          if (flag == 1)
          {
            // Loop over the nodes
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Derivatives w.r.t. solid positions will be handled by FDing
              /// Lagrange residual depends on interpolated pressure
              local_unknown = nodal_local_eqn(l2, P_index_interface);


              if (local_unknown >= 0)
              {
                // Jack - I have added a term here
                // This is for interpolated pressure in residual equation

                if (Q_zero == 1)
                {
                  // Full equations if analysis in Q_inv
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * psif(l) * W * J;

                  // Full equations is analysis in Q
                  //	      jacobian(local_eqn,local_unknown) +=
                  // Q*psif[l2]*psif(l)*W*J;
                }

                if (Q_zero == 0)
                {
                  // Q_inv = 0
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * psif(l) * W * J;

                  // Q = 0
                  // jacobian(local_eqn,local_unknown) +=
                  // 0;
                }
              }

              /// Lagrange residual also depends on smoothed derivatives through
              /// curvature term
              for (unsigned i2 = 0; i2 < 2; i2++)
              {
                local_unknown = projected_tangent_deriv_local_eqn(l2, i2);
                if (local_unknown >= 0)
                {
                  // Jack - I have taken out a factor of 1/Q
                  // This bit relates to d_kappa /d tau.
                  if (Q_zero == 1)
                  {
                    // Full equations if analysis in Q_inv
                    jacobian(local_eqn, local_unknown) +=
                      psif(l) * W * J * psif(l2) * interpolated_n[i2] /
                      (3.0 * alpha * alpha * Q);

                    // Full equations if analysis in Q
                    //  jacobian(local_eqn,local_unknown) +=
                    // psif(l)*W*J*psif(l2)*interpolated_n[i2]/(3.0*alpha*alpha);
                  }
                  if (Q_zero == 0)
                  {
                    // Q_inv = 0
                    jacobian(local_eqn, local_unknown) += 0;

                    // Q = 0
                    // jacobian(local_eqn,local_unknown) +=
                    // psif(l)*W*J*psif(l2)*interpolated_n[i2]/(3.0*alpha*alpha);
                  }
                }
              }
            }
          }
        }

        // Eqn for Lagrange multiplier (tangential lagrange constraint)
        //-------------------------------------------------------------
        local_eqn = tangential_lagrange_local_eqn(l);
        if (local_eqn >= 0)
        {
          residuals[local_eqn] += tangential_speed * psif(l) * W * J * sign;
          //                std::cout << "Tangential lagrange is not pinned" <<
          //                std::endl;

          if (flag == 1)
          {
            /// No nodal data in this residual!
          }
        }

        // Contribution to arbitrary solution measure
        //---------------------------------------------
        if (measure_eqn >= 0)
        {
          //                residuals[measure_eqn] += W * J * kappa*kappa;
          residuals[measure_eqn] += 0;
          if (flag == 1)
          {
            /// Currently a null measure.
          }
        }

        /// JACK, add extra integral measures here
        // Contribution to integral measures
        unsigned integral_values =
          external_data_pt(integral_measure_data_number)->nvalue();

        for (unsigned i_measure = 0; i_measure < integral_values; i_measure++)
        {
          double Q = 1 / Q_inv;
          double alpha = local_alpha;
          double ds = W * J * psif(l);
          double b = h;
          double pressure_inside_fluid = 0.5; // Q_zero*interpolated_p;
          double pressure_inside_bubble =
            0.5; // Q_zero*p_bubble + 1.0/(3.0*alpha*b*Q);//Jack

          if ((Q_zero != 1) && (Q_zero != 0))
          {
            std::cout << " " << std::endl;
            std::cout << "Oh oh, this ain't gonna work - need to set_Q_zero to "
                         "either 0 or 1 and you haven't"
                      << std::endl;
            exit(1);
          }

          if (Q_zero == 1)
          {
            // Full equations if analysis in Q_inv
            pressure_inside_fluid = interpolated_p;
            pressure_inside_bubble = p_bubble + 1.0 / (3.0 * alpha * b * Q);

            // Full equations if analysis in Q
            //		    pressure_inside_fluid = Q*interpolated_p;
            // pressure_inside_bubble = Q*p_bubble + 1.0/(3.0*alpha*b);
          }
          if (Q_zero == 0)
          {
            // Set Q_inv = 0
            pressure_inside_fluid = interpolated_p;
            pressure_inside_bubble = p_bubble;

            // Set Q = 0
            // pressure_inside_fluid = 0;
            // pressure_inside_bubble = 1.0/(3.0*alpha*b);
          }

          //                std::cout << pressure_inside_bubble << std::endl;
          double x = interpolated_x[0];
          double y = interpolated_x[1];
          /// Components of a normal vector which points out of the bubble
          double n_x = -interpolated_n[0];
          double n_y = -interpolated_n[1];


          local_eqn =
            external_local_eqn(integral_measure_data_number, i_measure);
          if (local_eqn >= 0)
          {
            if (i_measure == 0)
            {
              /// The perimeter of the bubble
              residuals[local_eqn] += 1 * ds;
            }
            else if (i_measure == 1)
            {
              /// The contribution to surface area from the vertical sides
              residuals[local_eqn] += b * ds;
            }
            else if (i_measure == 2)
            {
              /// The projected area of the bubble can be written as an integral
              /// over the interface by using the divergence theorem (or
              /// equivalently Green's theorem)
              residuals[local_eqn] += x * n_x * ds;
              /// residuals[local_eqn] += y*n_y*ds;
            }
            else if (i_measure == 3)
            {
              /// The volume of the bubble is also a line integral
              residuals[local_eqn] += b * x * n_x * ds;
            }
            else if (i_measure == 4)
            {
              /// We can construct the total surface area: 2 projected areas,
              /// and 1 vertical sides
              residuals[local_eqn] += (2 * x * n_x + b / alpha) * ds;
            }
            else if (i_measure == 5)
            {
              residuals[local_eqn] += y * b * ds;
            }
            else if (i_measure == 6)
            {
              /// Fluid pressure
              residuals[local_eqn] += pressure_inside_fluid * ds;
            }
            else if (i_measure == 7)
            {
              /// Bubble pressure
              residuals[local_eqn] += pressure_inside_bubble * ds;
            }
            else if (i_measure == 8)
            {
              /// Should integrate to zero
              residuals[local_eqn] += n_x * ds;
            }
            else if (i_measure == 9)
            {
              /// Should integrate to zero
              residuals[local_eqn] += n_y * ds;
            }
            else if (i_measure == 10)
            {
              /// This should integrate to 2*pi
              residuals[local_eqn] += kappa * ds;
            }
            else if (i_measure == 11)
            {
              /// This should integrate to zero!

              // If analysis in Q_inv
              residuals[local_eqn] +=
                (pressure_inside_fluid - pressure_inside_bubble +
                 1.0 / (3.0 * alpha * Q) * (1.0 / b + kappa / alpha)) *
                ds;

              // If analysis in Q
              // residuals[local_eqn] +=
              // (pressure_inside_fluid-pressure_inside_bubble+1.0/(3.0*alpha)*(1.0/b+kappa/alpha))*ds;
            }
            else if (i_measure == 12)
            {
              /// JACK - THIS SHOULD INTEGRATE TO THE AREA OF THE BUBBLE ALSO -
              /// Actually using Green's Theorem Area = (1/2)*integral_{bounding
              /// curve} xdx + ydy. Therefore Area = integral_{bounding curve}
              /// xdx = integral_{bounding_curve}ydy. Not int_{curve} y dy =
              /// int_{curve} y*dy/ds ds = int_{curve} y*n_y ds as required.
              residuals[local_eqn] += y * n_y * ds;
            }
          }
        }


        // Contribution to bulk equation (kinematic bc)
        //---------------------------------------------
        local_eqn = nodal_local_eqn(l, P_index_interface);

        if (local_eqn >= 0)
        {
          residuals[local_eqn] += h * St * normal_speed_time * psif(l) * W * J;
          residuals[local_eqn] += h * normal_speed_wall * psif(l) * W * J;

          /// These residuals depend only on solid position, and
          /// possibly external data. (but not nodal data)

          /// True mass matrix.
          if (flag == 2)
          {
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              /// Do x-components.
              local_unknown = this->position_local_eqn(l2, 0, 0);
              if (local_unknown >= 0)
              {
                mass_matrix(local_eqn, local_unknown) +=
                  h * St * psif(l) * W * J * interpolated_n[0] * psif(l2);
              }

              /// Do y-components
              local_unknown = this->position_local_eqn(l2, 0, 1);
              if (local_unknown >= 0)
              {
                mass_matrix(local_eqn, local_unknown) +=
                  h * St * psif(l) * W * J * interpolated_n[1] * psif(l2);
              }
            }
          }
        }

        // Lagrange multiplier contributions to pseudo_solid equations
        //------------------------------------------------------------
        for (unsigned i = 0; i < 2; i++)
        {
          local_eqn = this->position_local_eqn(l, 0, i);
          //                local_eqn = this->position_local_eqn(l,i);

          if (local_eqn >= 0)
          {
            // Add in a "Lagrange multiplier"
            residuals[local_eqn] +=
              -interpolated_lagrange * interpolated_n[i] * psif[l] * W * J;

            residuals[local_eqn] += interpolated_tangential_lagrange *
                                    interpolated_tangent[i] * psif[l] * W * J *
                                    sign;

            ////                    Do the Jacobian calculation
            if (flag == 1)
            {
              // Loop over the nodes
              for (unsigned l2 = 0; l2 < n_node; l2++)
              {
                // Derivatives w.r.t. solid positions will be handled by FDing
                // That leaves the "lagrange multipliers" only
                local_unknown = lagrange_local_eqn(l2);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    -psif[l2] * interpolated_n[i] * psif[l] * W * J;
                }
                local_unknown = tangential_lagrange_local_eqn(l2);
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    psif[l2] * interpolated_tangent[i] * psif[l] * W * J;
                }
              }
            } // End of Jacobian calculation
          }
        }

        /// Constructing the right multiplier R so that R is full rank
        /// and M*R is symmetric (and more sparse than M).
        if (flag == 3)
        {
          for (unsigned l2 = 0; l2 < n_node; l2++)
          {
            /// Do x-components. This is the pressure column.
            local_eqn = this->nodal_local_eqn(l, P_index_interface);
            local_unknown = this->position_local_eqn(l2, 0, 0);
            if (local_unknown >= 0 && local_eqn >= 0)
            {
              mass_matrix(local_unknown, local_eqn) +=
                h * St * psif(l) * W * J * interpolated_n[0] * psif(l2);
              //                        std::cout << "Adding " <<
              //                        h*St*psif(l)*W*J*interpolated_n[0]*psif(l2)
              //                        << "to " << local_unknown << " " <<
              //                        local_eqn << std::endl;
              ////                        mass_matrix(local_unknown,local_eqn)
              ////                        += 1;
            }


            /// Do y-components. This is the pressure column.
            local_eqn = this->nodal_local_eqn(l, P_index_interface);
            local_unknown = this->position_local_eqn(l2, 0, 1);
            if (local_unknown >= 0 && local_eqn >= 0)
            {
              mass_matrix(local_unknown, local_eqn) +=
                h * St * psif(l) * W * J * interpolated_n[1] * psif(l2);

              //                        mass_matrix(local_unknown,local_eqn)+=
              //                        1;
            }

            /// Do x-components. This is the pressure column.
            local_eqn = this->position_local_eqn(l, 0, 0);
            local_unknown = this->position_local_eqn(l2, 0, 0);
            if (local_unknown >= 0 && local_eqn >= 0)
            {
              mass_matrix(local_unknown, local_eqn) +=
                h * St * psif(l) * W * J * interpolated_n[1] * psif(l2);

              //                        mass_matrix(local_unknown,local_eqn)+=
              //                        1;
            }


            /// Do y-components. This is the pressure column.
            local_eqn = this->position_local_eqn(l, 0, 0);
            local_unknown = this->position_local_eqn(l2, 0, 1);
            if (local_unknown >= 0 && local_eqn >= 0)
            {
              mass_matrix(local_unknown, local_eqn) -=
                h * St * psif(l) * W * J * interpolated_n[0] * psif(l2);

              //                        mass_matrix(local_unknown,local_eqn)+=
              //                        1;
            }

            local_eqn = this->position_local_eqn(l, 0, 1);
            local_unknown = this->nodal_local_eqn(l2, P_index_interface);
            if (local_unknown >= 0 && local_eqn >= 0 && l == l2)
            {
              //                        std::cout << "Adding 1 to " << local_eqn
              //                        << " " << local_unknown << std::endl;
              mass_matrix(local_unknown, local_eqn) = 1;

              //                        mass_matrix(local_unknown,local_eqn)+=
              //                        1;
            }
          }
        }

      } // End of loop over shape functions
    } // End of loop over integration points


    if (flag == 3)
    {
      /// We also need to fill in some 'I' entries.

      for (unsigned l = 0; l < n_node; l++)
      {
        /// Lagrange multiplier
        local_eqn = this->nodal_local_eqn(l, 1);
        if (local_eqn >= 0)
        {
          mass_matrix(local_eqn, local_eqn) += 7;
        }
        /// Xdd equation
        local_eqn = this->nodal_local_eqn(l, 2);
        if (local_eqn >= 0)
        {
          mass_matrix(local_eqn, local_eqn) += 5;
        }
        /// Xdd equation
        local_eqn = this->nodal_local_eqn(l, 3);
        if (local_eqn >= 0)
        {
          mass_matrix(local_eqn, local_eqn) += 4.5;
        }
      }
    }


    if (flag == 1)
    {
      if (*fill_solid_by_fd_pt == true)
      {
        this->fill_in_jacobian_from_solid_position_by_fd(jacobian);
      }

      if (*fill_solid_by_fd_pt == false)
      {
        /////
        unsigned long a, b = 2, c = 3;
        a = this->ndof();
        RankThreeTensor<double> d_residual_dnodal(a, b, c, 0);

        get_dresidual_dnodal_coordinates(d_residual_dnodal);
        //
        for (unsigned n = 0; n < 3; n++)
        {
          for (unsigned d = 0; d < 2; d++)
          {
            local_eqn = this->position_local_eqn(n, 0, d);
            if (local_eqn >= 0)
            {
              for (unsigned u = 0; u < a; u++)
              {
                jacobian(u, local_eqn) = d_residual_dnodal(u, d, n);
              }
            }
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


} // namespace oomph

template<class ELEMENT>
void HeleShawInterfaceElement<ELEMENT>::get_dresidual_dnodal_coordinates(
  RankThreeTensor<double>& dresidual_dnodal_coordinates)
{
  const unsigned n_node = this->nnode();

  /// Each node has two coordinates.
  unsigned DIM = 2;
  // Set up memory for the shape and test functions
  Shape psif(n_node);
  DShape dpsifds(n_node, 1);
  DenseMatrix<double> dJ_dX(DIM, n_node);

  // Determine the number of integration points
  const unsigned n_intpt = this->integral_pt()->nweight();
  //
  // Integer to store the local equation number
  int local_eqn = 0;

  // Loop over the integration points
  for (unsigned ipt = 0; ipt < n_intpt; ipt++)
  {
    // Get the integral weight
    double W = this->integral_pt()->weight(ipt);

    this->dshape_local_at_knot(ipt, psif, dpsifds);

    // This is all we need! Jacobian is calculated from tangent vector.

    // We need to sort out:
    Vector<double> interpolated_x(DIM, 0.0);
    Vector<double> interpolated_t1(DIM, 0.0);

    double interpolated_Xdd = 0;
    double interpolated_Ydd = 0;
    double interpolated_p = 0;
    double interpolated_lagrange = 0;


    // We want to know how the two components of x
    // at each of the three nodes'
    // affect t1[0], t1[1] and y
    // at this particular ipt.

    DenseMatrix<double> d_t1_0_dX(DIM, n_node, 0.0);
    DenseMatrix<double> d_t1_1_dX(DIM, n_node, 0.0);
    DenseMatrix<double> dy_dX(DIM, n_node, 0.0);
    DenseMatrix<double> dn0_dX(DIM, n_node, 0.0);
    DenseMatrix<double> dn1_dX(DIM, n_node, 0.0);
    // Derivative of Jacobian of mapping w.r.t. to nodal coords


    for (unsigned l = 0; l < n_node; l++) // Looping over shape function.
    {
      // Jack - nodal_position(l,0) -- x ordinate
      //     - nodal_position(l,1) -- y ordinate
      //     - nodal_value(l,0) -- pressure
      //     - nodal_value(l,1) -- lagrangian mult
      //     - nodal_value(l,2) -- Xdd - don't know - solid equations
      //     - nodal_value(l,3) -- Ydd - as above

      interpolated_x[0] += this->nodal_position(l, 0) * psif(l);
      interpolated_x[1] += this->nodal_position(l, 1) * psif(l);
      interpolated_t1[0] +=
        this->nodal_position(l, 0) * dpsifds(l, 0); // Tangent vector
      interpolated_t1[1] +=
        this->nodal_position(l, 1) * dpsifds(l, 0); // Tangent vector
      interpolated_p += this->nodal_value(l, 0) * psif(l);
      interpolated_lagrange += this->nodal_value(l, 1) * psif(l);
      interpolated_Xdd += this->nodal_value(l, 2) * psif(l);
      interpolated_Ydd += this->nodal_value(l, 3) * psif(l);

      dy_dX(0, l) = 0.0;
      dy_dX(1, l) = psif(l);

      // These look like d/dX(t1) = d/ds(psif)

      d_t1_0_dX(0, l) += dpsifds(l, 0);
      d_t1_0_dX(1, l) += 0.0;

      d_t1_1_dX(0, l) += 0.0;
      d_t1_1_dX(1, l) += dpsifds(l, 0);
    }

    double tlength = interpolated_t1[0] * interpolated_t1[0] +
                     interpolated_t1[1] * interpolated_t1[1];

    /// Fill in dJ_dX.
    double J = std::sqrt(tlength);
    for (unsigned p = 0; p < DIM; p++)
    {
      for (unsigned q = 0; q < n_node; q++)
      {
        dJ_dX(p, q) = 0;
        dJ_dX(p, q) += interpolated_t1[0] * d_t1_0_dX(p, q) / J;
        dJ_dX(p, q) += interpolated_t1[1] * d_t1_1_dX(p, q) / J;
      }
    }

    /// FD to find dn_dX.
    Vector<double> n_plus(2, 0.0), n_minus(2, 0.0);

    double fd_step = 1e-8;
    for (unsigned l = 0; l < n_node; l++)
    {
      for (unsigned p = 0; p < DIM; p++)
      {
        double* const value_pt = &(node_pt(l)->x_gen(0, p));
        const double old_var = *value_pt;

        *value_pt += fd_step;
        outer_unit_normal(ipt, n_plus);
        *value_pt -= 2 * fd_step;
        outer_unit_normal(ipt, n_minus);

        dn0_dX(p, l) = (n_plus[0] - n_minus[0]) / (2 * fd_step);
        dn1_dX(p, l) = (n_plus[1] - n_minus[1]) / (2 * fd_step);

        *value_pt = old_var;
      }
    }

    // Now calculate the unit normal vector

    Vector<double> interpolated_n(2);
    outer_unit_normal(ipt, interpolated_n);

    Vector<double> U_wall(2, 0.0);
    get_wall_velocity(interpolated_x, U_wall);

    Vector<double> tau(2, 0.0);
    // Loop over the shape functions
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over directional components
      for (unsigned i = 0; i < 2; i++)
      {
        // Smoothed tangent vector
        tau[i] += projected_tangent_deriv(l, i) * psif(l);
      }
    }

    double kappa = 0.0;
    double normal_speed_wall = 0.0;
    for (unsigned k = 0; k < 2; k++)
    {
      kappa += tau[k] * interpolated_n[k];
      normal_speed_wall += U_wall[k] * interpolated_n[k];
    }

    double h = 0.0;
    double dhdt = 0.0;
    Vector<double> dhdx(2, 0.0);
    Vector<double> d_dhdt_dx(2, 0.0);

    // Jack - no longer seg faults but need to fill in the actual details in the
    // custom_hele_shaw etc. file
    upper_wall_flux_fct_pt()(interpolated_x, h, dhdt, dhdx, d_dhdt_dx);

    double b;
    double db_dy;
    b = h;
    db_dy = dhdx[1];

    // Get the value of the bubble pressure. This need not be constant.
    double bubble_pressure = get_p_bubble(interpolated_x);

    // This gives us dp/dy
    Vector<double> grad_pb(2, 0.0);
    {
      double y_init = interpolated_x[1];
      double pb_plus = 0, pb_minus = 0;
      interpolated_x[1] += 1e-8;
      pb_plus = get_p_bubble(interpolated_x);
      interpolated_x[1] -= 2e-8;
      pb_minus = get_p_bubble(interpolated_x);
      grad_pb[1] = (pb_plus - pb_minus) / 2e-8;
      interpolated_x[1] = y_init;
    }

    double alpha = this->alpha();

    double Q = 1.0 / this->q_inv();

    // Assemble d res_{local_eqn} / d X_{pq}
    // -------------------------------------

    // We're still at a fixed ipt.
    // Loop over the test functions
    for (unsigned l = 0; l < n_node; l++)
    {
      /// Loop over coordinate directions
      for (unsigned p = 0; p < DIM; p++)
      {
        /// Loop over nodes
        for (unsigned q = 0; q < n_node; q++)
        {
          // We have to deal with how x and y nodal components affect:

          // To determine Lagrange multiplier.
          // Depends on J, even if Ca_inv = 0.


          // JACK
          // We have the following equations
          //
          //
          //

          // First scope is basically differentiating (2) w.r.t X - the nodal
          // coordinates(which come as part of the solution)
          //


          //  residuals[local_eqn] +=(interpolated_p-p_bubble +
          //  kappa/(3.0*alpha*alpha*Q))*psif(l)*W*J;               (2)

          local_eqn = this->nodal_local_eqn(l, 1);

          if (local_eqn >= 0)
          {
            // Firstly, terms due to Jacobian varying (but no curvature changes)
            // TICK
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              W * psif(l) *
              (interpolated_p - bubble_pressure +
               kappa / (3.0 * alpha * alpha * Q)) *
              dJ_dX(p, q); // This comes from differentiating J w.r.t X in (2)

            // Terms due to pb varying- JACK - why no extra term?
            // This comes from differentiating p_bubble w.r.t X - grad_pb =
            // dp/dX = dp/dy*dy/dX + dp/dx*dx/dX;

            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              -W * psif(l) * J * grad_pb[1] * dy_dX(p, q);

            // This looks like W*psi*J*Ca_inv_k* (Xdd \dot dn/dX)

            // These two equations are for when we do d/dX(kappa)
            //
            //    d/dX(kappa) = d/dX[ d/ds(t) \dot n ]
            //                = d/dX[d/ds(t)] \dot n + d/ds(t) \dot d/dX(n)
            //
            //
            // This is the second term from above
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              W * psif(l) * J * 1.0 / (3.0 * alpha * alpha * Q) *
              (interpolated_Xdd * dn0_dX(p, q) +
               interpolated_Ydd *
                 dn1_dX(p, q)); // What is this? What is interpolated Xdd?
            // Terms due to curvature - not needed if curvature is in bubble
            // pressure.
            // This looks like -W*psi*Ca_inv_k*(Xdd \dot dt/dX)

            // This is the first term from above
            // dresidual_dnodal_coordinates(local_eqn,p,q) +=
            //-W*psif(l)*1.0/(3.0*alpha*alpha*Q)
            //*(interpolated_Xdd*d_t1_1_dX(p,q) -
            //interpolated_Ydd*d_t1_0_dX(p,q))/J; //What is this? - probably the
            //curavture term

            //			  dresidual_dnodal_coordinates(local_eqn,p,q) = 100;
          }

          // TODO: Time derivatives.

          //  residuals[local_eqn] +=
          //  h*sum(U_wall[k]*interpolated_n[k])*psif(l)*W*J;            (4)

          // These contributions probably come from (4)

          local_eqn = this->nodal_local_eqn(l, 0);

          if (local_eqn >= 0)
          {
            // looks like U_b*psif*W*J*b*d/dX(t) - JACK - do we need this term?

            // dresidual_dnodal_coordinates(local_eqn,p,q) +=
            // U_wall[0]*psif(l)*W*d_t1_1_dX(p,q)*b;

            // looks like U_b*psif*W*J*b*d/dX(n) - tick don't need the U_wall[1]
            // component because it is zero
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              U_wall[0] * psif(l) * W * b * dn0_dX(p, q) * J;

            // looks like U_b*n*psif*W*J*d/dX(b) - tick

            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              U_wall[0] * interpolated_n[0] * psif(l) * W * db_dy *
              dy_dX(p, q) * J;

            // looks like U_b*n*psif*W*b*dJ/dX - tick

            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              U_wall[0] * interpolated_n[0] * psif(l) * W * b * dJ_dX(p, q);
          }

          local_eqn = this->nodal_local_eqn(l, 2);

          // Only tau is not based on a solid position...apparantly
          //  residuals[local_eqn]
          //  +=(tau[i]*psif(l)*J+interpolated_tangent[i]*dpsifds(l,0))*W; (1)

          if (local_eqn >= 0)
          {
            // TICK
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              W * dJ_dX(p, q) * interpolated_Xdd * psif(l);

            // looks like W*d/dX(t)*d/ds(psif)/J -> Why /J ?
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              W * d_t1_0_dX(p, q) * dpsifds(l, 0) / J;

            // looks like -W*t*d/dX(J)*d/ds(psif)/J^2 -> Why /J*J ?
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              -W * interpolated_t1[0] * dJ_dX(p, q) * dpsifds(l, 0) / (J * J);
          }

          local_eqn = this->nodal_local_eqn(l, 3);

          if (local_eqn >= 0)
          {
            // As last scope but for Ydd

            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              W * dJ_dX(p, q) * interpolated_Ydd * psif(l);

            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              W * d_t1_1_dX(p, q) * dpsifds(l, 0) / J;

            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              -W * interpolated_t1[1] * dJ_dX(p, q) * dpsifds(l, 0) / (J * J);
          }


          ///////////////////////////////////////////////////////////////////////////////////
          // These next two are for the solid equations
          // TICK
          //  residuals[local_eqn]
          //  +=-interpolated_lagrange*interpolated_n[i]*psif[l]*W*J; (5)

          // Jack - where are the contributuons to the equation in (6)?

          local_eqn = this->position_local_eqn(l, 0, 0);

          if (local_eqn >= 0)
          {
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              -interpolated_lagrange * interpolated_n[0] * psif(l) * W *
              dJ_dX(p, q); // Differentiating the J term in (5)

            dresidual_dnodal_coordinates(
              local_eqn, p, q) += // Differentiating the interpolated n in (5)
              -interpolated_lagrange * dn0_dX(p, q) * W * psif(l) * J;
          }

          local_eqn = this->position_local_eqn(l, 0, 1);

          if (local_eqn >= 0)
          {
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              -interpolated_lagrange * interpolated_n[1] * psif(l) * W *
              dJ_dX(p, q); // As above

            dresidual_dnodal_coordinates(local_eqn, p, q) +=
              -interpolated_lagrange * dn1_dX(p, q) * W * psif(l) * J;
          }
          //////////////////////////////////////////////////////////////////////////////////
        }
      }
    }
  } // End of loop over integration points

} /// END DRESIDUAL
