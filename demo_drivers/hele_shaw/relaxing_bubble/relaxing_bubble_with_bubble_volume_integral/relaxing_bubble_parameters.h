#ifndef OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER
#define OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER

#include "generic.h"

namespace relaxing_bubble
{
  double ca_inv;
  double st;
  double alpha;
  /// Pseudo-solid Poisson ratio
  double nu = 0;
  double* bubble_pressure_pt = 0;
  double* target_bubble_volume_pt = 0;

  ConstitutiveLaw* constitutive_law_pt = 0;

  void depth_fct(const Vector<double>& x, double& f)
  {
    f = 1.0;
  }

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    depth_fct(x, b);
    dbdt = 0.0;
  }

  void wall_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    U_frame[0] = 0.0;
    U_frame[1] = 0.0;
  }

  void bubble_pressure_fct(const Vector<double>& x, double& pressure)
  {
    pressure = (*bubble_pressure_pt);
  }
} // namespace relaxing_bubble

#endif
