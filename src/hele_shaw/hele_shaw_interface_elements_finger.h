#ifndef OOMPH_HELE_SHAW_INTERFACE_ELEMENTS_FINGER_HEADER
#define OOMPH_HELE_SHAW_INTERFACE_ELEMENTS_FINGER_HEADER


namespace oomph
{
  //=======================================================================
  /// 1D Free surface elements for Hele Shaw problems with pseudo-elastic
  /// mesh motion
  //======================================================================
  template<class ELEMENT>
  class HeleShawInterfaceElementFinger : public virtual FaceGeometry<ELEMENT>,
                                         public virtual FaceElement
  {
  public:
    /// \short Function pointer to function which provides bubble pressure as a
    /// function of vector x. This should usually be a constant - should default
    /// to zero. Any perturbations to the system can be supplied through a
    /// spatially varying bubble pressure function.
    typedef void (*BubblePressureFctPt)(const Vector<double>& x,
                                        double& p_bubble);

    /// \short Function pointer to function which provides wall speed as a
    /// function of x. This allows us to solve for bubble motion in a moving
    /// frame. A constant wall speed does not affect the mass conservation
    /// equations, but does feature in the kinematic equation for interface
    /// motion.
    typedef void (*WallSpeedFctPt)(const Vector<double>& x,
                                   Vector<double>& U_wall);

    /// Access function: Pointer to bubble pressure function
    BubblePressureFctPt& bubble_pressure_fct_pt()
    {
      return Bubble_pressure_fct_pt;
    }

    /// Access function: Pointer to bubble pressure function. Const version
    BubblePressureFctPt bubble_pressure_fct_pt() const
    {
      return Bubble_pressure_fct_pt;
    }

    /// Access function: Pointer to moving frame speed function
    WallSpeedFctPt& wall_speed_fct_pt()
    {
      return Wall_speed_fct_pt;
    }

    /// Access function: Pointer to moving frame speed function. Const version
    WallSpeedFctPt wall_speed_fct_pt() const
    {
      return Wall_speed_fct_pt;
    }

    /// Function to set all Lagrange multipliers in the element to zero
    void set_lagrange_multipliers_to_zero()
    {
      unsigned n_node = this->nnode();
      for (unsigned inod = 0; inod < n_node; inod++)
      {
        lagrange(inod) = 0.0;
      }
    }

  private:
    /// Flag to indicate whether we apply thin-film effects or not
    bool Apply_thin_film_effects;

    /// Thin film effect contribution to the kinematic BC
    double thin_film_effect_kinematic_bc(double& Ca)
    {
      if (Apply_thin_film_effects)
      {
        return thin_film_homotopy() * pow(Ca, 2.0 / 3.0) /
               (0.76 + 2.16 * pow(Ca, 2.0 / 3.0));
      }
      else
      {
        return 0.0;
      }
    }

    /// Thin film effect contribution to the dynamic BC
    double thin_film_effect_dynamic_bc(double& Ca)
    {
      if (Apply_thin_film_effects)
      {
        return 1.0 +
               thin_film_homotopy() *
                 (pow(Ca, 2.0 / 3.0) / (0.26 + 1.48 * pow(Ca, 2.0 / 3.0)) +
                  1.59 * Ca);
      }
      else
      {
        return 1.0;
      }
    }

    /// Pointer to function that specifies the bubble pressure function
    BubblePressureFctPt Bubble_pressure_fct_pt;

    /// Pointer to function that specifies the moving frame speed function
    WallSpeedFctPt Wall_speed_fct_pt;

    /// \short Pointer to the aspect_ratio ratio: reference gap width /
    /// in-plane lengthscale
    double* Aspect_ratio_pt;

    /// \short Pointer to the finger tip node
    Node* Finger_tip_node_pt;

    Vector<double>* Local_coordinate_finger_tip_pt;

    //   GeneralisedElement* Element_containing_tip_node_pt;
    ELEMENT** Element_containing_tip_node_pt_pt;

    /// \short Pointer to membrane height at bubble tip
    double* Membrane_height_at_tip_pt;

    double* Ha_pt;

    // \short Pointer to the distance to the opposite point at the interface
    Vector<double>* Vector_of_distance_to_opposite_point_pt;

    Vector<double>* Vector_of_ipt_index_of_opposite_points_pt;

    Vector<HeleShawInterfaceElementFinger<ELEMENT>*>*
      Vector_of_pointers_to_opposite_elements_pt;

    /// \short Pointer to the fjord index
    unsigned* Fjord_index_pt;

    /// Pointer to the inverse Capillary number
    double* Ca_inv_pt;

    /// Pointer to the Strouhal number
    double* St_pt;

    /// Pointer to a homotopy parameter to control thin film effects
    double* Thin_film_homotopy_pt;

    /// Pointer to a homotopy parameter to control the thickness of the film
    /// effects
    double* Thickness_homotopy_pt;

    /// Default value for physical constants
    static double Default_Physical_Constant_Value;

    /// \short ID of additional unknowns introduced by this face element
    /// (smoothed components of derivative of tangent vector, and
    /// Lagrange multiplier)
    unsigned Id;

    /// Index of the nodal value at which the pressure is stored
    unsigned P_index_interface;

    /// Equation number of the equation for the Lagrange multiplier
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
    /// for the derivative of the i-th component of the tangent vector at node j
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

    void get_curvature_at_the_interface(
      Vector<Vector<double>>& vector_of_curvature);

  public:
    /// \short Constructor, pass a pointer to the bulk element and the face
    /// index of the bulk element to which the element is to be attached to.
    /// The optional identifier can be used to distinguish the additional nodal
    /// values created by this element (in order:
    /// 0: Lagrange multiplier;
    /// 1: x-component of smoothed derivative of tangent vector;
    /// 2: y-component of smoothed derivative of tangent vector)
    /// from those created by other FaceElements.
    HeleShawInterfaceElementFinger(FiniteElement* const& element_pt,
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

      // Set the Strouhal number to the default value
      St_pt = &Default_Physical_Constant_Value;

      // Find the index at which the pressure stored
      // from the bulk element
      ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(element_pt);
      this->P_index_interface = bulk_element_pt->p_index_hele_shaw();

      // Read out the number of nodes on the face
      unsigned n_node_face = this->nnode();

      // Set the additional data values in the face
      // There are three additional values at each node -- the
      // Lagrange multiplier and the two components of the smoothed
      // derivative of the tangent vector
      Vector<unsigned> additional_data_values(n_node_face);
      for (unsigned i = 0; i < n_node_face; i++)
      {
        additional_data_values[i] = 3;
      }

      // Now add storage and set the map containing
      // the position of the first entry of this face element's
      // additional values.
      add_additional_values(additional_data_values, id);

      // By default we don't apply thin film effects
      Apply_thin_film_effects = false;
      Thin_film_homotopy_pt = 0;
      Thickness_homotopy_pt = 0;
    }


    // rmodif
    /// Return the contribution to the integral of the square of the deviation
    /// from the mean radius
    void compute_deviation_from_mean_radius(const double& mean_radius,
                                            double& deviation_contribution)
    {
      // Initialise
      deviation_contribution = 0.0;

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local coordinate
      Vector<double> s(1);

      // Storage for the global coordinate
      Vector<double> x(2);

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
        Vector<double> interpolated_tangent(2, 0.0);
        Vector<double> interpolated_x(2, 0.0);

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            // Spatial bits
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_tangent[i] +=
              this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Calculate the length of the tangent Vector
        double tlength = interpolated_tangent[0] * interpolated_tangent[0] +
                         interpolated_tangent[1] * interpolated_tangent[1];

        // Set the Jacobian of the line element
        double J = sqrt(tlength);

        // Add contribution
        deviation_contribution += W * J *
                                  (sqrt(interpolated_x[0] * interpolated_x[0] +
                                        interpolated_x[1] * interpolated_x[1]) -
                                   mean_radius) *
                                  (sqrt(interpolated_x[0] * interpolated_x[0] +
                                        interpolated_x[1] * interpolated_x[1]) -
                                   mean_radius);
      }
    }

    // rmodif
    /// Return the contribution to the mean radius
    void compute_mean_radius(double& arc_length_contribution,
                             double& mean_radius_nod_contribution)
    {
      // Initialise
      mean_radius_nod_contribution = 0.0;
      arc_length_contribution = 0.0;

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local coordinate
      Vector<double> s(1);

      // Storage for the global coordinate
      Vector<double> x(2);

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
        Vector<double> interpolated_tangent(2, 0.0);
        Vector<double> interpolated_x(2, 0.0);

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            // Spatial bits
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_tangent[i] +=
              this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Calculate the length of the tangent Vector
        double tlength = interpolated_tangent[0] * interpolated_tangent[0] +
                         interpolated_tangent[1] * interpolated_tangent[1];

        // Set the Jacobian of the line element
        double J = sqrt(tlength);

        // Add contribution
        arc_length_contribution += W * J;
        mean_radius_nod_contribution +=
          W * J *
          sqrt(interpolated_x[0] * interpolated_x[0] +
               interpolated_x[1] * interpolated_x[1]);
      }
    }


    ///////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////// ////Joao
    ///=====>>>>> added function
    void compute_mean_y(double& arc_length_contribution, double& mean_y_nod)
    {
      // Initialise
      mean_y_nod = 0.0;
      arc_length_contribution = 0.0;

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local coordinate
      Vector<double> s(1);

      // Storage for the global coordinate
      Vector<double> x(2);

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
        Vector<double> interpolated_tangent(2, 0.0);
        Vector<double> interpolated_x(2, 0.0);

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            // Spatial bits
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_tangent[i] +=
              this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Calculate the length of the tangent Vector
        double tlength = interpolated_tangent[0] * interpolated_tangent[0] +
                         interpolated_tangent[1] * interpolated_tangent[1];

        // Set the Jacobian of the line element
        double J = sqrt(tlength);

        // Add contribution
        arc_length_contribution += W * J; ////// W*J
        mean_y_nod += W * J * interpolated_x[1]; ///// W*J*interpolated_x[1]
      }
      mean_y_nod /= arc_length_contribution;
    }

    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////

    /// Return the contribution to the flow into the film
    double compute_flow_into_thin_films()
    {
      // Initialise
      double flow_into_thin_films = 0.0;

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();

      // Get the value of the inverse Capillary number
      double Ca_inv = ca_inv();

      // Get the values of the membrane height at the bubble tip
      //     double h_tip = membrane_height_at_tip();
      double h_tip = get_height_at_tip();

      // Get the value of the Strouhal numer
      double St = st();

      // Storage for the local cooridinate
      Vector<double> s(1);
      Vector<double> s_bulk(2);

      // Get pointer to bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

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
        Vector<double> interpolated_x(2, 0.0);
        Vector<double> interpolated_dx_dt(2, 0.0);
        Vector<double> interpolated_tangent(2, 0.0);

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            // Spatial bits
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_dx_dt[i] += this->dnodal_position_dt(l, i) * psif(l);
            interpolated_tangent[i] +=
              this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Get h from bulk!
        get_local_coordinate_in_bulk(s, s_bulk);

        // Get the velocity of the plates
        Vector<double> U_wall(2, 0.0);
        get_wall_velocity(interpolated_x, U_wall);

        // Get the pressure gradient
        Vector<double> pressure_gradient(2, 0.0);
        bulk_el_pt->get_pressure_gradient(s_bulk, pressure_gradient);

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

        // Speed of the bubble boundary
        double normal_speed_interface = 0.0;

        // Horizontal speed of the plates
        double normal_speed_wall = 0.0;

        double normal_pressure_gradient = 0.0;

        for (unsigned k = 0; k < 2; k++)
        {
          normal_speed_interface += interpolated_dx_dt[k] * interpolated_n[k];
          normal_speed_wall += U_wall[k] * interpolated_n[k];
          normal_pressure_gradient += pressure_gradient[k] * interpolated_n[k];
        }

        // double interface_velocity =
        // sqrt(interpolated_dx_dt[0]*interpolated_dx_dt[0] +
        //      interpolated_dx_dt[1]*interpolated_dx_dt[1]) +
        // sqrt(U_wall[0]*U_wall[0] + U_wall[1]*U_wall[1]);

        // Modulus of the frame velocity ////Joao
        double U_frame_modulus =
          sqrt(U_wall[0] * U_wall[0] + U_wall[1] * U_wall[1]);

        // obacht
        double Ca_local = U_frame_modulus / Ca_inv;
        //       double Ca_local = 1.0/Ca_inv; // * interface_velocity;


        double h = 0.0;
        double dhdt = 0.0;
        Vector<double> dhdx(2, 0.0);

        //---------------------------------------------------------
        // NOTE: This is dangerous/wrong. We're really supposed to
        // pass in dpsi/dx (2D!) which we don't have/need here.
        // However it's only used for the computation of the ALE bits
        // in dh/dt and dh/dt isn't used here! Phew...
        //---------------------------------------------------------
        // Number of nodes in bulk element
        unsigned n_nod_bulk = bulk_el_pt->nnode();
        DShape dpsidx_bulk(n_nod_bulk, 2);
        Shape psi_bulk(n_nod_bulk);
        bulk_el_pt->shape(s_bulk, psi_bulk);
        bulk_el_pt->get_upper_wall_data(
          s_bulk, interpolated_x, psi_bulk, dpsidx_bulk, h, dhdt, dhdx);

        double effective_h =
          thickness_homotopy() * h_tip + (1.0 - thickness_homotopy()) * h;

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Contribution to the flow into the thin films
          //---------------------------------------------
          //         flow_into_thin_films -= h*St*normal_speed_interface*
          //          thin_film_effect_kinematic_bc(Ca_local)*psif(l)*W*J;
          flow_into_thin_films -= effective_h * St * normal_speed_interface *
                                  thin_film_effect_kinematic_bc(Ca_local) *
                                  psif(l) * W * J;


          // obacht, we assume zero pressure gradient in the films
          //         flow_into_thin_films -= pow(h,3.0)*
          //          pow(thin_film_effect_kinematic_bc(Ca_local),3.0)*
          //          normal_pressure_gradient*psif(l)*W*J;
          flow_into_thin_films -=
            pow(effective_h, 3.0) *
            pow(thin_film_effect_kinematic_bc(Ca_local), 3.0) *
            normal_pressure_gradient * psif(l) * W * J;

          //         flow_into_thin_films -= h*normal_speed_wall*
          //          thin_film_effect_kinematic_bc(Ca_local)*psif(l)*W*J;
          flow_into_thin_films -= effective_h * normal_speed_wall *
                                  thin_film_effect_kinematic_bc(Ca_local) *
                                  psif(l) * W * J;
        }
      }

      return flow_into_thin_films;
    }


    ///\short Compute the element's residual vector and the Jacobian matrix.
    /// Jacobian is computed by finite-differencing.
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Add the contribution to the residuals
      this->fill_in_contribution_to_residuals(residuals);

      // Allocate storage for the full residuals (residuals of entire element)
      unsigned n_dof = this->ndof();
      Vector<double> full_residuals(n_dof);

      // Get the residuals for the entire element
      this->get_residuals(full_residuals);

      // Get the solid entries in the jacobian using finite differences
      // oomph_info<<"In face element before"<<std::endl;
      this->fill_in_jacobian_from_solid_position_by_fd(full_residuals,
                                                       jacobian);
      // oomph_info<<"In face element after"<<std::endl;

      // There could be internal data
      //(finite-difference the lot by default)
      this->fill_in_jacobian_from_internal_by_fd(
        full_residuals, jacobian, true);

      // There could also be external data
      //(finite-difference the lot by default)
      // oomph_info<<"Contribution from BulkElement"<<std::endl;
      this->fill_in_jacobian_from_external_by_fd(
        full_residuals, jacobian, true);

      // There could also be nodal data
      this->fill_in_jacobian_from_nodal_by_fd(full_residuals, jacobian);
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

    void fill_in_contribution_to_mass_matrix(Vector<double>& residuals,
                                             DenseMatrix<double>& mass_matrix)
    {
      // Add the residual contributions
      oomph_info << "Filli in contribution to mass matrix in "
                    "HeleShawInterfaceElementFinger ...\n";
      fill_in_generic_residual_contribution_hele_shaw_interface(
        residuals, GeneralisedElement::Dummy_matrix, mass_matrix, 2);
      oomph_info << "Filling in contribution to mass matrix in "
                    "HeleShawInterfaceElementFinger done!\n";
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Add the residual contributions
      oomph_info << "Filli in contribution to jacobian and mass matrix in "
                    "HeleShawInterfaceElementFinger ...\n";
      fill_in_generic_residual_contribution_hele_shaw_interface(
        residuals, GeneralisedElement::Dummy_matrix, mass_matrix, 2);

      // Allocate storage for the full residuals (residuals of entire element)
      unsigned n_dof = this->ndof();
      Vector<double> full_residuals(n_dof);

      // Get the residuals for the entire element
      this->get_residuals(full_residuals);

      // Get the solid entries in the jacobian using finite differences
      // oomph_info<<"In face element before"<<std::endl;
      this->fill_in_jacobian_from_solid_position_by_fd(full_residuals,
                                                       jacobian);
      // oomph_info<<"In face element after"<<std::endl;

      // There could be internal data
      //(finite-difference the lot by default)
      this->fill_in_jacobian_from_internal_by_fd(
        full_residuals, jacobian, true);

      // There could also be external data
      //(finite-difference the lot by default)
      // oomph_info<<"Contribution from BulkElement"<<std::endl;
      this->fill_in_jacobian_from_external_by_fd(
        full_residuals, jacobian, true);

      // There could also be nodal data
      this->fill_in_jacobian_from_nodal_by_fd(full_residuals, jacobian);
      oomph_info << "Filling in contribution tojacobian and mass matrix in "
                    "HeleShawInterfaceElementFinger done!\n";
    }

    void get_curvature(Vector<Vector<double>>& vector_of_curvature)
    {
      get_curvature_at_the_interface(vector_of_curvature);
    }

    inline void update_in_external_fd(const unsigned& i)
    {
      // Get the current bubble pressure
      double bubble_pressure = this->external_data_pt(0)->value(0);
      // Loop over all nodes
      unsigned n_node = this->nnode();
      for (unsigned inod = 0; inod < n_node; inod++)
      {
        Node* nod_pt = this->node_pt(inod);
        // Only update when the node is in the interior of the bubble where the
        // pressure dof is pinned
        if (nod_pt->is_pinned(4))
        {
          // oomph_info << "Updating pressure!"<<std::endl;
          nod_pt->set_value(4, bubble_pressure);
        }
      }
    }

    inline void update_before_external_fd()
    {
      const unsigned i = 0;
      update_in_external_fd(i);
    }

    inline void update_before_nodal_fd()
    {
      const unsigned i = 0;
      update_in_external_fd(i);
    }

    inline void update_before_solid_position_fd()
    {
      const unsigned i = 0;
      update_in_external_fd(i);
    }

    /// function to enable thin film effects
    void enable_thin_film_effects()
    {
      Apply_thin_film_effects = true;
    }

    /// function to enable thin film effects
    void disable_thin_film_effects()
    {
      Apply_thin_film_effects = false;
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
    const double& ca_inv() const
    {
      return *Ca_inv_pt;
    }

    /// Pointer to the inverse Capillary number
    double*& ca_inv_pt()
    {
      return Ca_inv_pt;
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

    /// Aspect ratio: Reference gap width / in-plane lengthscale
    double aspect_ratio() const
    {
#ifdef PARANOID
      if (Aspect_ratio_pt == 0)
      {
        throw OomphLibError(
          "Aspect_ratio_pt has not been set yet for HeleShaw elements",
          "HeleShawInterfaceElementFinger::aspect_ratio()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Aspect_ratio_pt;
    }

    /// \short Pointer to aspect ratio
    double*& aspect_ratio_pt()
    {
      return Aspect_ratio_pt;
    }

    /// \short Pointer to aspect ratio. Const version.
    double* aspect_ratio_pt() const
    {
      return Aspect_ratio_pt;
    }


    /// Finger tip node
    Node finger_tip_node() const
    {
#ifdef PARANOID
      if (Finger_tip_node_pt == 0)
      {
        throw OomphLibError(
          "Finger_tip_node_pt has not been set yet for HeleShaw elements",
          "HeleShawInterfaceElementFinger::finger_tip_node()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Finger_tip_node_pt;
    }

    /// \short Pointer to aspect ratio
    Node*& finger_tip_node_pt()
    {
      return Finger_tip_node_pt;
    }

    /// \short Pointer to aspect ratio. Const version.
    Node* finger_tip_node_pt() const
    {
      return Finger_tip_node_pt;
    }


    /// Local coordinate of finger tip
    Vector<double> local_coordinate_finger_tip() const
    {
#ifdef PARANOID
      if (Local_coordinate_finger_tip_pt == 0)
      {
        throw OomphLibError(
          "Local_coordinate_finger_tip_pt has not been set yet for HeleShaw "
          "elements",
          "HeleShawInterfaceElementFinger::local_coordinate_finger_tip()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Local_coordinate_finger_tip_pt;
    }

    /// \short Pointer to local coordinate of finger tip
    Vector<double>*& local_coordinate_finger_tip_pt()
    {
      return Local_coordinate_finger_tip_pt;
    }

    /// \short Pointer to local coordinate of finger tip. Const version.
    Vector<double>* local_coordinate_finger_tip_pt() const
    {
      return Local_coordinate_finger_tip_pt;
    }


    /// Local coordinate of finger tip
    //    GeneralisedElement element_containing_tip_node() const
    ELEMENT element_containing_tip_node() const
    {
#ifdef PARANOID
      if (Element_containing_tip_node_pt_pt == 0)
      {
        throw OomphLibError(
          "Element_containing_tip_node_pt_pt has not been set yet for HeleShaw "
          "elements",
          "HeleShawInterfaceElementFinger::element_containing_tip_node()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return **Element_containing_tip_node_pt_pt;
    }

    /// \short Pointer to local coordinate of finger tip
    //    GeneralisedElement*& element_containing_tip_node_pt()
    ELEMENT**& element_containing_tip_node_pt_pt()
    {
      return Element_containing_tip_node_pt_pt;
    }

    /// \short Pointer to local coordinate of finger tip. Const version.
    //    GeneralisedElement* element_containing_tip_node_pt() const
    ELEMENT** element_containing_tip_node_pt_pt() const
    {
      return Element_containing_tip_node_pt_pt;
    }

    ELEMENT* element_containing_tip_node_pt() const
    {
      return *Element_containing_tip_node_pt_pt;
    }


    double get_height_at_tip()
    {
      double height_at_tip = 0.0;

      height_at_tip =
        1.0 + (element_containing_tip_node_pt()->interpolated_w_fvk(
                local_coordinate_finger_tip(), 0)) /
                aspect_ratio();

      return height_at_tip;
    }
    //    double get_height_at_tip()
    //     {
    //      Vector<double> x_auxiliary(2);
    //      Vector<double> s_auxiliary(2);
    //
    //      double height_at_tip = 0.0;
    //
    //      if(finger_tip_node_pt() != 0)
    //       {
    //        x_auxiliary[0]=finger_tip_node_pt()->x(0);
    //        x_auxiliary[1]=finger_tip_node_pt()->x(1);
    //       }
    //      else
    //       {
    //        x_auxiliary[0]=0.0;
    //        x_auxiliary[1]=0.0;
    //       }
    //
    //      // Sub-geomobject (FE) that contains that point
    //      GeomObject* sub_geom_object_pt=0;
    //
    //      // Find it
    //      Mesh_as_geom_object_pt->locate_zeta(x_auxiliary, sub_geom_object_pt,
    //      s_auxiliary);
    //
    //      if(sub_geom_object_pt!=0)
    //       {
    //        // Cast to our element
    //        ELEMENT* el_pt=dynamic_cast<ELEMENT*>(sub_geom_object_pt);
    //        if (el_pt!=0)
    //         {
    //          height_at_tip =
    //            1.0 +
    //            (el_pt->interpolated_w_fvk(s_auxiliary,0))/aspect_ratio();
    //         }
    //        else
    //         {
    //          throw OomphLibError(
    //          "el_pt was not set",
    //          "HeleShawInterfaceElementFinger::get_height_at_tip()",
    //          OOMPH_EXCEPTION_LOCATION);
    //          abort();
    //         }
    //       }
    //      else
    //       {
    //        throw OomphLibError(
    //        "sub_geom_object_pt was not set",
    //        "HeleShawInterfaceElementFinger::get_height_at_tip()",
    //        OOMPH_EXCEPTION_LOCATION);
    //        abort();
    //       }
    //
    //      return height_at_tip;
    //     }

    /// Membrane height at bubble tip
    double membrane_height_at_tip() const
    {
#ifdef PARANOID
      if (Membrane_height_at_tip_pt == 0)
      {
        throw OomphLibError(
          "Membrane_height_at_tip_pt has not been set yet for HeleShaw "
          "elements",
          "HeleShawInterfaceElementFinger::membrane_height_at_tip()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Membrane_height_at_tip_pt;
    }


    /// \short Pointer to membrane height at bubble tip
    double*& membrane_height_at_tip_pt()
    {
      return Membrane_height_at_tip_pt;
    }

    /// \short Pointer to membrane height at bubble tip. Const version.
    double* membrane_height_at_tip_pt() const
    {
      return Membrane_height_at_tip_pt;
    }


    Vector<double> vector_of_ipt_index_of_opposite_points() const
    {
#ifdef PARANOID
      if (Vector_of_ipt_index_of_opposite_points_pt == 0)
      {
        throw OomphLibError(
          "Vector_of_ipt_index_of_opposite_points_pt has not been set yet for "
          "HeleShaw elements",
          "HeleShawInterfaceElementFinger::vector_of_ipt_index_of_opposite_"
          "points()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Vector_of_ipt_index_of_opposite_points_pt;
    }

    Vector<double>*& vector_of_ipt_index_of_opposite_points_pt()
    {
      return Vector_of_ipt_index_of_opposite_points_pt;
    }

    Vector<double>* vector_of_ipt_index_of_opposite_points_pt() const
    {
      return Vector_of_ipt_index_of_opposite_points_pt;
    }

    // vector_of_pointers_to_opposite_elements_pt()
    // Vector_of_pointers_to_opposite_elements_pt
    // Vector< HeleShawInterfaceElementFinger<ELEMENT>* >*

    Vector<HeleShawInterfaceElementFinger<ELEMENT>*> vector_of_pointers_to_opposite_elements()
      const
    {
#ifdef PARANOID
      if (Vector_of_pointers_to_opposite_elements_pt == 0)
      {
        throw OomphLibError(
          "Vector_of_pointers_to_opposite_elements_pt has not been set yet for "
          "HeleShaw elements",
          "HeleShawInterfaceElementFinger::vector_of_pointers_to_opposite_"
          "elements()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Vector_of_pointers_to_opposite_elements_pt;
    }

    Vector<HeleShawInterfaceElementFinger<ELEMENT>*>*& vector_of_pointers_to_opposite_elements_pt()
    {
      return Vector_of_pointers_to_opposite_elements_pt;
    }

    Vector<HeleShawInterfaceElementFinger<ELEMENT>*>* vector_of_pointers_to_opposite_elements_pt()
      const
    {
      return Vector_of_pointers_to_opposite_elements_pt;
    }


    double ha() const
    {
#ifdef PARANOID
      if (Ha_pt == 0)
      {
        throw OomphLibError("Ha_pt has not been set yet for HeleShaw elements",
                            "HeleShawInterfaceElementFinger::ha()",
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Ha_pt;
    }

    double*& ha_pt()
    {
      return Ha_pt;
    }

    double* ha_pt() const
    {
      return Ha_pt;
    }


    Vector<double> vector_of_distance_to_opposite_point() const
    {
#ifdef PARANOID
      if (Vector_of_distance_to_opposite_point_pt == 0)
      {
        throw OomphLibError(
          "Vector_of_distance_to_opposite_point_pt has not been set yet for "
          "HeleShaw elements",
          "HeleShawInterfaceElementFinger::vector_of_distance_to_opposite_"
          "point()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Vector_of_distance_to_opposite_point_pt;
    }

    Vector<double>*& vector_of_distance_to_opposite_point_pt()
    {
      return Vector_of_distance_to_opposite_point_pt;
    }

    Vector<double>* vector_of_distance_to_opposite_point_pt() const
    {
      return Vector_of_distance_to_opposite_point_pt;
    }


    /// Fjord index
    unsigned fjord_index() const
    {
#ifdef PARANOID
      if (Fjord_index_pt == 0)
      {
        throw OomphLibError(
          "Fjord_index_pt has not been set yet for HeleShaw elements",
          "HeleShawInterfaceElementFinger::fjord_index()",
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Fjord_index_pt;
    }

    /// \short Pointer to fjord index
    unsigned*& fjord_index_pt()
    {
      return Fjord_index_pt;
    }

    /// \short Pointer to fjord index. Const version.
    unsigned* fjord_index_pt() const
    {
      return Fjord_index_pt;
    }

    /// Thin film homotopy parameter
    double thin_film_homotopy() const
    {
      if (Thin_film_homotopy_pt == 0)
      {
        return 1.0;
      }
      else
      {
        return *Thin_film_homotopy_pt;
      }
    }

    /// \short Pointer to thin film homotopy parameter
    double*& thin_film_homotopy_pt()
    {
      return Thin_film_homotopy_pt;
    }

    /// \short Pointer to thin film homotopy parameter. Const version.
    double* thin_film_homotopy_pt() const
    {
      return Thin_film_homotopy_pt;
    }


    /// Thin thickness homotopy parameter
    double thickness_homotopy() const
    {
      if (Thickness_homotopy_pt == 0)
      {
        return 1.0;
      }
      else
      {
        return *Thickness_homotopy_pt;
      }
    }

    /// \short Pointer to thickness homotopy parameter.
    double*& thickness_homotopy_pt()
    {
      return Thickness_homotopy_pt;
    }

    /// \short Pointer to thickness homotopy parameter. Const version.
    double* thickness_homotopy_pt() const
    {
      return Thickness_homotopy_pt;
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

    /// Return the value of the horizontal speed of the plates
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

      // Get the value of the inverse Capillary number
      double Ca_inv = ca_inv();

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

        double interpolated_lagrange = 0.0;

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Lagrange multiplier
          interpolated_lagrange += this->lagrange(l) * psif(l);

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
          sigma_kappa += sigma_local * Ca_inv * tau[k] * interpolated_n[k];
        }

        outfile << interpolated_x[0] << " " << interpolated_x[1] << " "
                << tau[0] << " " << tau[1] << " " << interpolated_dx_dt[0]
                << " " << interpolated_dx_dt[1] << " " << sigma_kappa << " "
                << interpolated_lagrange << " " << std::endl;
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
      // hierher fill this in when we've finalised what lives in output
    }
  };


  //============================================================
  /// Default value for physical constant (static)
  //============================================================
  template<class ELEMENT>
  double
    HeleShawInterfaceElementFinger<ELEMENT>::Default_Physical_Constant_Value =
      1.0;


  //=======================================================================
  /// Calculate the residual contribution (kinematic and dynamic BC and
  /// Lagrange multiplier contribution from pseudo-elastic node updates
  /// Flag here is not the same as the "standard" oomph-lib convention.
  /// Instead we have 0 for just residuals, 1 for residuals + jacobian
  /// residuals + mass matrix will be calculated for flag value of 2
  //=======================================================================
  template<class ELEMENT>
  void HeleShawInterfaceElementFinger<ELEMENT>::
    fill_in_generic_residual_contribution_hele_shaw_interface(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, 1);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Get the value of the inverse Capillary number
    double Ca_inv = ca_inv();

    // Get the value of the membrane height at the bubble tip
    //   double h_tip = membrane_height_at_tip();
    double h_tip = get_height_at_tip();


    double Ha = ha();

    Vector<double> vector_opposite_distances =
      vector_of_distance_to_opposite_point();
    Vector<double> vector_ipt_index_of_opposite_points =
      vector_of_ipt_index_of_opposite_points();

    Vector<HeleShawInterfaceElementFinger<ELEMENT>*>
      vector_of_pointers_to_opposite_elements_in_the_fjord =
        vector_of_pointers_to_opposite_elements();

    // Get the value of the Strouhal numer
    double St = st();

    // Integers to store the local equation numbers
    int local_eqn = 0;
    int local_unknown;

    // Storage for the local cooridinate
    Vector<double> s(1);

    // Get pointer to bulk element
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

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
      Vector<double> interpolated_tangent(2, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      Vector<double> interpolated_dx_dt(2, 0.0);
      Vector<double> tau(2, 0.0);

      // Loop over the shape functions
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_lagrange += lagrange(l) * psif[l];
        interpolated_p += node_pt(l)->value(P_index_interface) * psif[l];

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

      // Get h from bulk!
      Vector<double> s_bulk(2);
      get_local_coordinate_in_bulk(s, s_bulk);

      // Get the value of the bubble pressure
      double p_bubble = get_p_bubble(interpolated_x);

      // Get the velocity of the plates
      Vector<double> U_wall(2, 0.0);
      get_wall_velocity(interpolated_x, U_wall);

      // Get the pressure gradient
      Vector<double> pressure_gradient(2, 0.0);
      bulk_el_pt->get_pressure_gradient(s_bulk, pressure_gradient);

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

      // Also get the (possibly variable) surface tension
      // double sigma_local = this->sigma(s);

      // Assemble the surface tension and normal speed terms
      // double sigma_kappa = 0.0;

      // Curvature of the bubble
      double kappa = 0.0;

      // Speed of the bubble boundary
      double normal_speed_interface = 0.0;

      // Horizontal speed of the plates
      double normal_speed_wall = 0.0;

      double normal_pressure_gradient = 0.0;

      for (unsigned k = 0; k < 2; k++)
      {
        // sigma_kappa += sigma_local*Ca_inv*tau[k]*interpolated_n[k];
        kappa += tau[k] * interpolated_n[k];
        normal_speed_interface += interpolated_dx_dt[k] * interpolated_n[k];
        normal_speed_wall += U_wall[k] * interpolated_n[k];
        normal_pressure_gradient += pressure_gradient[k] * interpolated_n[k];
      }

      // double interface_velocity =
      // sqrt(interpolated_dx_dt[0]*interpolated_dx_dt[0] +
      //      interpolated_dx_dt[1]*interpolated_dx_dt[1]) +
      // sqrt(U_wall[0]*U_wall[0] + U_wall[1]*U_wall[1]);

      // Modulus of the frame velocity ////////Joao
      double U_frame_modulus =
        sqrt(U_wall[0] * U_wall[0] +
             U_wall[1] * U_wall[1]); /// Joao => Multiply Ca by U_frame_modulus,
                                     /// so the speed of the bubble is used to
                                     /// calculate the thin film effect

      // obacht
      double Ca_local =
        U_frame_modulus /
        Ca_inv; /// Joao => Multiply Ca by U_frame_modulus, so the speed of the
                /// bubble is used to calculate the thin film effect
      //     double Ca_local = 1.0/Ca_inv; // * interface_velocity;


      double h = 0.0;
      double dhdt = 0.0;
      Vector<double> dhdx(2, 0.0);

      //---------------------------------------------------------
      // NOTE: This is dangerous/wrong. We're really supposed to
      // pass in dpsi/dx (2D!) which we don't have/need here.
      // However it's only used for the computation of the ALE bits
      // in dh/dt and dh/dt isn't used here! Phew...
      //---------------------------------------------------------
      // Number of nodes in bulk element
      unsigned n_nod_bulk = bulk_el_pt->nnode();
      DShape dpsidx_bulk(n_nod_bulk, 2);
      Shape psi_bulk(n_nod_bulk);
      bulk_el_pt->shape(s_bulk, psi_bulk);
      bulk_el_pt->get_upper_wall_data(
        s_bulk, interpolated_x, psi_bulk, dpsidx_bulk, h, dhdt, dhdx);

      double effective_h =
        thickness_homotopy() * h_tip + (1.0 - thickness_homotopy()) * h;

      // Non-dim gap-width
      double local_aspect_ratio = aspect_ratio();

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
                jacobian(local_eqn, local_unknown) +=
                  psif[l2] * psif(l) * W * J;
              }
            }
          } // End of Jacobian calculation
        }


        // Eqn for Lagrange multiplier (dynamic free surface condition)
        //-------------------------------------------------------------
        local_eqn = lagrange_local_eqn(l);
        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            (interpolated_p - p_bubble +
             local_aspect_ratio * Ca_inv *
               (2.0 * thin_film_effect_dynamic_bc(Ca_local) / h +
                local_aspect_ratio * kappa) /
               12.0) *
            psif(l) * W * J;


          /////////////////////////////////////////////////////
          //// Testing disjoining pressure
          // vector_of_pointers_to_opposite_elements_in_the_fjord
          // vector_ipt_index_of_opposite_points
          /////////////////////////////////////////////////////
          if (Ha != 0.0)
          {
            if (vector_ipt_index_of_opposite_points[ipt] > -1.0)
            {
              double opposite_ipt =
                unsigned(vector_ipt_index_of_opposite_points[ipt]);

              Vector<Vector<double>> opposite_element_information;
              vector_of_pointers_to_opposite_elements_in_the_fjord[ipt]
                ->get_curvature_at_the_interface(opposite_element_information);

              Vector<double> opposite_interpolated_x(2);
              opposite_interpolated_x[0] =
                opposite_element_information[opposite_ipt][0];
              opposite_interpolated_x[1] =
                opposite_element_information[opposite_ipt][1];

              //             double fjord_length_scale =
              //             h_tip/local_aspect_ratio;
              double fjord_length_scale = 1.0 * 1.0e-2;
              double fjord_thickness =
                sqrt(pow(interpolated_x[0] - opposite_interpolated_x[0], 2.0) +
                     pow(interpolated_x[1] - opposite_interpolated_x[1], 2.0));

              //             double disjoining_p =
              //               (Ha/(pow(fjord_thickness,3.0)))*(1.0 -
              //               pow(fjord_length_scale/fjord_thickness,3.0));

              //             double disjoining_p =
              //             -Ha/(6.0*3.14159*pow(fjord_thickness,3.0));
              double disjoining_p =
                -Ha * pow(fjord_length_scale / fjord_thickness, 3.0);
              residuals[local_eqn] += -(disjoining_p)*psif(l) * W * J;
            }
          }


          /////////////////////////////////////////////////////


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
                // This is for interpolated pressure in residual equation
                jacobian(local_eqn, local_unknown) +=
                  psif[l2] * psif(l) * W * J;
              }

              /// Lagrange residual also depends on smoothed derivatives through
              /// curvature term
              for (unsigned i2 = 0; i2 < 2; i2++)
              {
                local_unknown = projected_tangent_deriv_local_eqn(l2, i2);
                if (local_unknown >= 0)
                {
                  // This bit relates to d_kappa /d tau.
                  jacobian(local_eqn, local_unknown) +=
                    local_aspect_ratio * Ca_inv * local_aspect_ratio / 12.0 *
                    psif(l) * W * J * psif(l2) * interpolated_n[i2];
                }
              }
            }
          }
        }

        // Contribution to bulk equation (kinematic bc)
        //---------------------------------------------
        local_eqn = nodal_local_eqn(l, P_index_interface);
        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            h * St * normal_speed_interface * psif(l) * W * J;

          // Moving frame terms
          residuals[local_eqn] += h * normal_speed_wall * psif(l) * W * J;

          // Thin film contributions
          //         residuals[local_eqn] -= h*St*normal_speed_interface*
          //          thin_film_effect_kinematic_bc(Ca_local)*psif(l)*W*J;
          residuals[local_eqn] -= effective_h * St * normal_speed_interface *
                                  thin_film_effect_kinematic_bc(Ca_local) *
                                  psif(l) * W * J;

          // obacht, we assume zero pressure gradient in the film
          //         residuals[local_eqn] -=
          //         pow(h,3.0)*normal_pressure_gradient*
          //          pow(thin_film_effect_kinematic_bc(Ca_local),3.0)*psif(l)*W*J;
          residuals[local_eqn] -=
            pow(effective_h, 3.0) * normal_pressure_gradient *
            pow(thin_film_effect_kinematic_bc(Ca_local), 3.0) * psif(l) * W * J;

          //         residuals[local_eqn] -= h*normal_speed_wall*
          //          thin_film_effect_kinematic_bc(Ca_local)*psif(l)*W*J;
          residuals[local_eqn] -= effective_h * normal_speed_wall *
                                  thin_film_effect_kinematic_bc(Ca_local) *
                                  psif(l) * W * J;


          if (flag == 2)
          {
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              /// Do x-components.
              local_unknown = this->position_local_eqn(l2, 0, 0);
              if (local_unknown >= 0)
              {
                mass_matrix(local_eqn, local_unknown) +=
                  h * St * interpolated_n[0] * psif(l) * psif(l2) * W * J;

                // Thin film contributions
                //               mass_matrix(local_eqn, local_unknown)
                //                -=
                //                h*St*thin_film_effect_kinematic_bc(Ca_local)*
                //                interpolated_n[0]*psif(l)*psif(l2)*W*J;
                mass_matrix(local_eqn, local_unknown) -=
                  effective_h * St * thin_film_effect_kinematic_bc(Ca_local) *
                  interpolated_n[0] * psif(l) * psif(l2) * W * J;
              }

              /// Do y-components
              local_unknown = this->position_local_eqn(l2, 0, 1);
              if (local_unknown >= 0)
              {
                mass_matrix(local_eqn, local_unknown) +=
                  h * St * interpolated_n[1] * psif(l) * psif(l2) * W * J;

                // Thin film contributions
                //               mass_matrix(local_eqn, local_unknown)
                //                -=
                //                h*St*thin_film_effect_kinematic_bc(Ca_local)*
                //                interpolated_n[1]*psif(l)*psif(l2)*W*J;
                mass_matrix(local_eqn, local_unknown) -=
                  effective_h * St * thin_film_effect_kinematic_bc(Ca_local) *
                  interpolated_n[1] * psif(l) * psif(l2) * W * J;
              }
            }
          }

          /// These residuals depend only on solid position, and
          /// possibly external data.
        }

        // Lagrange multiplier contributions to pseudo_solid equations
        //------------------------------------------------------------
        for (unsigned i = 0; i < 2; i++)
        {
          local_eqn = this->position_local_eqn(l, 0, i);

          if (local_eqn >= 0)
          {
            // Add in a "Lagrange multiplier"
            residuals[local_eqn] +=
              -interpolated_lagrange * interpolated_n[i] * psif[l] * W * J;

            // Do the Jacobian calculation
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
              }
            } // End of Jacobian calculation
          }
        }
      } // End of loop over shape functions
    } // End of loop over integration points
  }

  template<class ELEMENT>
  void HeleShawInterfaceElementFinger<ELEMENT>::get_curvature_at_the_interface(
    Vector<Vector<double>>& vector_of_curvature)
  {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, 1);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    vector_of_curvature.resize(n_intpt);

    // Integers to store the local equation numbers
    int local_eqn = 0;
    int local_unknown;

    // Storage for the local cooridinate
    Vector<double> s(1);

    // Get pointer to bulk element
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

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
      Vector<double> interpolated_tangent(2, 0.0);
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
          interpolated_tangent[i] += this->nodal_position(l, i) * dpsifds(l, 0);
        }
      }

      // Get h from bulk!
      //     Vector<double> s_bulk(2);
      //     get_local_coordinate_in_bulk(s,s_bulk);

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

      // Curvature of the bubble
      double kappa = 0.0;

      for (unsigned k = 0; k < 2; k++)
      {
        kappa += tau[k] * interpolated_n[k];
      }

      vector_of_curvature[ipt].resize(10);
      vector_of_curvature[ipt][0] = interpolated_x[0];
      vector_of_curvature[ipt][1] = interpolated_x[1];
      vector_of_curvature[ipt][2] = kappa;
      vector_of_curvature[ipt][3] = interpolated_n[0];
      vector_of_curvature[ipt][4] = interpolated_n[1];
      vector_of_curvature[ipt][5] = interpolated_tangent[0];
      vector_of_curvature[ipt][6] = interpolated_tangent[1];
      // Distance between the two sides of the interface in a fjord.
      // That is an arbitrary value >> 1, the real value is going
      // to be updated in the driver code (we need acess to all
      // surface elements to calculate it).
      vector_of_curvature[ipt][7] = 1000.0;
      // x position of opposit point.
      // That is an arbitrary value >> 1, the real value is going
      // to be updated in the driver code (we need acess to all
      // surface elements to calculate it).
      vector_of_curvature[ipt][8] = 1000.0;
      // y position of opposit point.
      // That is an arbitrary value >> 1, the real value is going
      // to be updated in the driver code (we need acess to all
      // surface elements to calculate it).
      vector_of_curvature[ipt][9] = 1000.0;

    } // End of loop over integration points

    /*
         Vector<double> node_coord(3);
         Vector<double> b_coord(1);
         unsorted_interface.resize(n_node);

         //loop over nodes
         for(unsigned inod=0; inod<n_node; inod++)
          {
           Node* nod_pt = this->node_pt(inod);

           nod_pt->get_coordinates_on_boundary(Inner_boundary,b_coord);
           node_coord[0] = b_coord[0];
           node_coord[1] = nod_pt->x(0);
           node_coord[2] = nod_pt->x(1);
           node_coord[3] = inod;

           unsorted_interface[inod].resize(4);
           unsorted_interface[inod] = node_coord;
          }
    */
  }


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


} // namespace oomph

#endif
