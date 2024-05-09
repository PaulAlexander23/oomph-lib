//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
namespace Problem_Parameter
{
  int global_Q_zero = 0;

  // Jack
  double hessian_tolerance;
  double tressian_tolerance;
  int hessian_order;
  int tressian_order;
  int eps_number = 0;
  int theta_number = 0;
  ofstream Trace_weakly_nonlinear;

  bool fill_in_solid_jacobian_interface_by_fd = false;

  bool hopf1;

  void channel_height_derivative(const Vector<double>& x, Vector<double>& dbdx)
  {
    // std::cout<<"I am here"<<std::endl;
    double y = x[1];
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

    double cosh_plus = pow(std::cosh(sharpness * (y + width)), 2);
    double cosh_minus = pow(std::cosh(sharpness * (y - width)), 2);

    dbdx[1] = -sharpness * height * 0.5 * (1.0 / cosh_plus - 1.0 / cosh_minus);
    dbdx[1] += -asymmetry;
    dbdx[0] = 0.0;
  }


  void upper_wall_flux_fct(const Vector<double>& x,
                           double& b,
                           double& dbdt,
                           Vector<double>& dbdx,
                           Vector<double>& d_dbdt_dx)
  {
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
      UpperWall_file << b << std::endl;
    }
    dbdt = 0.0;

    channel_height_derivative(x, dbdx);

    d_dbdt_dx[0] = 0.0;
    d_dbdt_dx[1] = 0.0;
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
    /// Constructor initialise error
    MyNewElement()
    {
      Error = 0.0;
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

        outfile << velocity[0] << " " << velocity[1] << " "
                << interpolated_p_hele_shaw(s) << " " << pressure_gradient[0]
                << " " << pressure_gradient[1] << " " << Error << " " << size()
                << " " << h << " "
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
