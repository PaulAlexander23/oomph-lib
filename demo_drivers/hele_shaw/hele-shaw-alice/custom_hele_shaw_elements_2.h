//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
namespace Problem_Parameter
{
  /// Doc info object
  DocInfo Doc_info;

  std::string Output_label = "test";

  /// Pseudo-solid Poisson ratio
  double Nu = 0.3;

  /// Initial radius of bubble
  double Radius = 0.7;

  bool Minisolve_state = false;
  bool ignore_height_effects_in_dynamic_bc = false;
  bool ignore_height_effects_in_viscous_terms = false;

  bool Set_steady_conditions_after_adapt = false;

  bool doc_every_step = false;

  double* global_G_pt = 0;
  double* global_G_upper_outflow_pt = 0;
  double* global_Q_inv_pt = 0;
  double* global_CoM_Y_pt = 0;
  double* global_Frame_speed_pt = 0;
  double* global_Obstacle_height_pt = 0;
  double* global_Obstacle_width_pt = 0;
  double* global_alpha_pt = 0;
  double* global_Asymmetry_pt = 0;
  double* global_St_pt = 0;
  double* global_delta_pt = 0;
  double* global_measure_pt = 0;
  double* global_drop_viscosity_pt = 0;

  unsigned branch_no = 0;

  bool Topology_change_needed = false;
  bool Valid_output_point = false;


  double Time_shift = 0;

  double Min_distance = 0;


  double x_center = 0.0;
  double y_center = 0.0;
  double e_vertical = 1.0; // The eccentricity of the initial circle

  unsigned N_Bubble = 1;
  Vector<double> Bubble_volume_vector(N_Bubble, 0);
  // Vector<double*> Bubble_pressure_vector_pt(N_Bubble);
  Vector<Vector<Vector<double>>> All_the_bubbles;
  Vector<Vector<Vector<double>>> Ordered_bubbles;
  Vector<double> Drop_viscosity_vector(N_Bubble, 0.1);

  Vector<double> vector_of_eigenvalues_rp(10, 2.0);
  Vector<double> vector_of_eigenvalues_ip(10, 2.0);

  bool Reload_from_vector = false;

  // Node* Tip_node_pt;
  bool ignore_bulk_effects = false;

  /// \short Volume of the bubble (negative because it's outside the
  /// fluid!)
  double Volume = -MathematicalConstants::Pi * Radius * Radius;
  double Centre_of_mass = 0;

  /// \short Scaling factor for inflow velocity (allows it to be switched off
  /// to do hydrostatics)
  double Inflow_veloc_magnitude = 0.0;

  /// \short Length of the channel
  double Length = 4.0;

  double Channel_start = -Length;
  double Channel_end = Length;

  bool Read_from_file = false;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw* Constitutive_law_pt = 0;

  /// Trace file
  ofstream Trace_file;

  ofstream M_file;

  /// \short File to document the norm of the solution (for validation
  /// purposes -- triangle doesn't give fully reproducible results so
  /// mesh generation/adaptation may generate slightly different numbers
  /// of elements on different machines!)
  ofstream Norm_file;

  void channel_height_function(const Vector<double>& x, double& b)
  {
    /// This function should have obstacle width and height, asymmetry and
    /// possibly sharpness.

    double height = 0.1;
    double width = 0.33;
    double asymmetry = 0;
    double sharpness = 40; /// Quite smooth for now.
    sharpness = 40; /// ALICE!
    if (global_Obstacle_height_pt != 0)
    {
      height = *global_Obstacle_height_pt;
    }
    if (global_Obstacle_width_pt != 0)
    {
      width = *global_Obstacle_width_pt;
    }
    if (global_Asymmetry_pt != 0)
    {
      asymmetry = *global_Asymmetry_pt;
    }
    double y = x[1];
    double tanh_plus = std::tanh(sharpness * (y + width));
    double tanh_minus = std::tanh(sharpness * (y - width));
    b = 1 - height * 0.5 * (tanh_plus - tanh_minus);
    b += -y * asymmetry;
  }

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
    }
    dbdt = 0;
  }

  void frame_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    if (global_Frame_speed_pt != 0)
    {
      U_frame[0] = *global_Frame_speed_pt;
      U_frame[1] = 0;
    }
    else
    {
      U_frame[0] = 0;
      U_frame[1] = 0;
    }
  }

  void drop_pressure_function(const Vector<double>& x, double& pressure)
  {
    // std::cout << i_bubble << std::endl;
    double height = 1.0;
    if (ignore_height_effects_in_dynamic_bc == false)
    {
      channel_height_function(x, height);
    }
    double alpha = *global_alpha_pt;
    double Q = 1 / (*global_Q_inv_pt);
    pressure = -1 / (3 * alpha * Q * height);
  }


  void bubble_pressure_gradient_function(const Vector<double>& x,
                                         Vector<double>& grad_pb)
  {
    /// ALICE_FLAG
    std::cout << "Evaluating bubble pressure gradient" << std::endl;
    std::cout << "Missing height dependence!" << std::endl;
  }


  void normal_flux_ahead_of_bubble(const Vector<double>& x, double& flux)
  {
    /// Ahead of the bubble we have the boundary condition dp/dx = -G.
    /// We need to supply the function b^3 n.grad p = -G b^3.
    double b;
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
    }
    double G;
    if (global_G_pt != 0)
    {
      G = *global_G_pt;
    }
    else
    {
      G = 1;
    }
    double p_x = -1 * G;
    flux = p_x * (b * b * b);
  }

  void normal_flux_behind_bubble(const Vector<double>& x, double& flux)
  {
    /// Ahead of the bubble we have the boundary condition dp/dx = -G.
    /// We need to supply the function b^3 n.grad p = -G b^3.
    double b;
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
    }
    double G;
    if (global_G_pt != 0)
    {
      G = *global_G_pt;
    }
    else
    {
      G = 1;
    }
    double p_x = G;
    flux = p_x * (b * b * b);
    //    std::cout << x[0] << " " << x[1] << " " << flux << std::endl;
  }

} // namespace Problem_Parameter


namespace oomph
{
  class MyNewElement
    : public virtual ProjectableHeleShawElement<
        PseudoSolidNodeUpdateElement<THeleShawElement<3>, TPVDElement<2, 3>>>
  {
  private:
    /// Storage for elemental error estimate -- used for post-processing
    double Error;

  public:
    double Viscosity;
    unsigned pressure_node_value;

    /// Constructor initialise error
    MyNewElement()
    {
      Error = 0.0;
      Viscosity = 1.0;
      HeleShawEquations::pressure_index = 1;
      /// We are hard coding to allow for up to 2 pressure values.
    }

    void set_pressure_node_value(int pressure_index)
    {
      HeleShawEquations::pressure_index = pressure_index;
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      /// The bulk element has no time derivative terms.
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Set error value for post-processing
    void set_error(const double& error)
    {
      Error = error;
    }

    void get_error(double& error)
    {
      error = Error;
    }

    /// Overload output function
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates and velocity
      Vector<double> s(2);
      Vector<double> x(2);

      Vector<double> velocity(2);
      Vector<double> pressure_gradient(2);

      double h, dhdt;
      Vector<double> dhdx(2), d_dhdt_dx(2);
      // Tecplot header info
      outfile << tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates and velocity at plot point
        get_s_plot(iplot, nplot, s);
        get_velocity(s, velocity);
        get_pressure_gradient(s, pressure_gradient);

        for (unsigned i = 0; i < 2; i++)
        {
          x[i] = interpolated_x(s, i);
          outfile << interpolated_x(s, i) << " ";
        }

        get_upper_wall_flux_data(0, x, h, dhdt, dhdx, d_dhdt_dx);

        outfile << interpolated_p_hele_shaw(s) << " "
                << interpolated_p_arbitrary(s, 0) << " "
                << interpolated_p_arbitrary(s, 1) << " " << Viscosity << " "
                << Error << " " << size() << " " << h << " "
                << double(pressure_index) << " " << velocity[0] << " "
                << velocity[1] << " "
                << "\n";
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      write_tecplot_zone_footer(outfile, nplot);
    }


    /// Overload output function
    void output(std::ostream& outfile)
    {
      unsigned nplot = 3;
      output(outfile, nplot);
    }
  };


  template<>
  class FaceGeometry<MyNewElement> : public virtual SolidTElement<1, 3>
  {
  public:
    FaceGeometry() : SolidTElement<1, 3>() {}
  };

  template<>
  class FaceGeometry<FaceGeometry<MyNewElement>>
    : public virtual SolidPointElement
  {
  public:
    FaceGeometry() : SolidPointElement() {}
  };


} // namespace oomph
