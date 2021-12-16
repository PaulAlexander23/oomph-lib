#ifndef OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER
#define OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER

#include "generic.h"

namespace moving_bubble
{
  /// Moving Bubble Problem Parameter --- Variables

  /// Dimensionless parameters
  double ca_inv;
  double st;
  double alpha;

  /// Total flux has be nondimensionalised to 1
  const double total_flux = 1.0;

  /// Bubble's initial condition
  double major_radius;
  double bubble_initial_centre_x = 1.0;
  double bubble_initial_centre_y;

  /// Channel perturbation parameters
  double perturbation_amplitude;
  double perturbation_rms_width;
  const double perturbation_centre_x = 1.5;
  const double perturbation_centre_y = 0.5;

  /// Volume pointers
  double target_bubble_volume = 0;
  double target_fluid_volume = 0;
  double total_volume = 0;

  /// X moment to constrain CoM
  double target_x_mom = 0;

  /// Pseudo-solid Poisson ratio pointer
  double nu;

  /// Reference frame speed
  double* global_frame_speed_pt = 0;
  double* global_frame_travel_pt = 0;

  /// Constitutive law pointer
  ConstitutiveLaw* constitutive_law_pt = 0;

  /// Bubble pressure pointer
  double* bubble_pressure_pt = 0;

  /// Flux b3 pointers
  double* inlet_b3_pt = 0;
  double* outlet_b3_pt = 0;

  /// Moving Bubble Problem Parameter --- Functions

  /// Channel depth without rate of change
  void channel_depth(const Vector<double>& x, double& b)
  {
    // Transform the domain as required i.e. from -1:1 to 0:1
    double local_x = x[0];
    double local_y = x[1];

    b = 1.0 -
        perturbation_amplitude *
          exp(-(local_x - (perturbation_centre_x - *global_frame_travel_pt)) *
                (local_x - (perturbation_centre_x - *global_frame_travel_pt)) /
                (2 * perturbation_rms_width * perturbation_rms_width) -
              (local_y - perturbation_centre_y) *
                (local_y - perturbation_centre_y) /
                (2 * perturbation_rms_width * perturbation_rms_width));
    //    b =
    //      1.0 - perturbation_amplitude / 2 *
    //              (cos(2 * MathematicalConstants::Pi * (local_x -
    //              perturbation_centre_x) / perturbation_rms_width) +
    //               cos(2 * MathematicalConstants::Pi * (local_y -
    //               perturbation_centre_y) / perturbation_rms_width));
  }

  /// Channel depth with rate of change
  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    channel_depth(x, b);
    dbdt = 0.0;
  }

  // Moving reference frame speed function
  void frame_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    U_frame[0] = (*global_frame_speed_pt);
    U_frame[1] = 0.0;
  }

  /// Bubble pressure function
  void bubble_pressure_fct(const Vector<double>& x, double& pressure)
  {
    pressure = (*bubble_pressure_pt);
  }

  /// Inlet flux boundary condition
  void get_inlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_b3 and total flux
    const double inlet_area = 1.0;
    double total_inlet_flux =
      total_flux - inlet_area * (*global_frame_speed_pt);
    double dpdx = total_inlet_flux / *inlet_b3_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b) / 12;
  }

  /// Outlet boundary condition
  void get_outlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_b3 and total flux
    const double outlet_area = 1.0;
    double total_outlet_flux =
      - total_flux + outlet_area * (*global_frame_speed_pt);
    double dpdx = total_outlet_flux / *outlet_b3_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b) / 12;
  }

  void print_parameters()
  {
    /// Dimensionless parameters
    cout << "ca_inv:" << ca_inv << endl;
    cout << "st:" << st << endl;
    cout << "alpha:" << alpha << endl;

    /// Total flux has be nondimensionalised to 1
    cout << "total_flux:" << total_flux << endl;

    /// Bubble's initial condition
    cout << "major_radius:" << major_radius << endl;
    cout << "bubble_initial_centre_x:" << bubble_initial_centre_x << endl;
    cout << "bubble_initial_centre_y:" << bubble_initial_centre_y << endl;

    /// Channel perturbation parameters
    cout << "perturbation_amplitude:" << perturbation_amplitude << endl;
    cout << "perturbation_rms_width:" << perturbation_rms_width << endl;
    cout << "perturbation_centre_x:" << perturbation_centre_x << endl;
    cout << "perturbation_centre_y:" << perturbation_centre_y << endl;


    /// Volume pointers
    cout << "target_bubble_volume:" << target_bubble_volume << endl;
    cout << "target_fluid_volume:" << target_fluid_volume << endl;
    cout << "total_volume:" << total_volume << endl;

    /// Pseudo-solid Poisson ratio pointer
    cout << "nu:" << nu << endl;
  }

  class ODEFunctor : public SolutionFunctorBase
  {
  public:
    /// Constructor
    ODEFunctor() {}

    /// Destructor
    virtual ~ODEFunctor() {}

    /// Exact or approximate solution. Used for initialisation and/or error
    /// checking
    Vector<double> operator()(const double& t, const Vector<double>& x) const
    {
      Vector<double> output(1);

      output[0] = (*global_frame_speed_pt) * t;

      return output;
    }

    /// Derivative function. Specifies the ODE that we are solving
    Vector<double> derivative(const double& t,
                              const Vector<double>& x,
                              const Vector<double>& u) const
    {
      Vector<double> output(1);

      output[0] = *global_frame_speed_pt;

      return output;
    }
  };
} // namespace moving_bubble

#endif
