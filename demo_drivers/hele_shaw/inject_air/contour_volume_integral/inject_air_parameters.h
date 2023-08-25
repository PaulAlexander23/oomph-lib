#ifndef OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER
#define OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER

#include "generic.h"

namespace inject_air
{
  /// Dimensionless parameters
  double ca_inv;
  double st;
  double alpha;

  /// Pseudo-solid Poisson ratio
  double nu;
  ConstitutiveLaw* constitutive_law_pt = 0;

  /// Bubble pressure pointer
  double* bubble_pressure_pt = 0;

  /// Bubble geometry
  const double circular_radius = 0.05;
  double major_radius;
  double* target_bubble_volume_pt = 0;
  double target_fluid_volume;
  double total_volume;

  const double initial_volume = MathematicalConstants::Pi * pow(0.05, 2.0);
  double injection_rate = 0.1;

  void channel_height_function(const Vector<double>& x, double& b)
  {
    b = 1.0;
  }

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    channel_height_function(x, b);
    dbdt = 0.0;
  }

  void wall_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    U_frame[0] = 0.0;
    U_frame[1] = 0.0;
  }

  void bubble_pressure_fct(const Vector<double>& x, double& pressure)
  {
    //cout << "bubble pressure:" << *bubble_pressure_pt << endl;
    pressure = *bubble_pressure_pt;
  }

  void volume_fct(const double& t, double& volume)
  {
    volume = - initial_volume - injection_rate * t;
  }

} // namespace inject_air

#endif
