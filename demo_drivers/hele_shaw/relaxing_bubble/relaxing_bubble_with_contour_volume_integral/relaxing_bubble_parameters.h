#ifndef OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER
#define OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER

#include "generic.h"

namespace relaxing_bubble
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
  const double circular_radius = 0.3;
  double major_radius;
  double target_bubble_volume;
  double target_fluid_volume;
  double total_volume;

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
} // namespace relaxing_bubble

#endif
