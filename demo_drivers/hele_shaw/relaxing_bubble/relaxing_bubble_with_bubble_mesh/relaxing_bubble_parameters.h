#ifndef OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER
#define OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER

#include "generic.h"

namespace relaxing_bubble
{
  double* q_inv_pt = 0;
  double* st_pt = 0;
  double* alpha_pt = 0;
  /// Pseudo-solid Poisson ratio
  double* nu_pt = 0;

  ConstitutiveLaw* constitutive_law_pt = 0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    b = 1.0;
    dbdt = 0.0;
  }

  void wall_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    U_frame[0] = 0.0;
    U_frame[1] = 0.0;
  }

  void bubble_pressure_fct(const Vector<double>& x, double& pressure)
  {
    double bubble_radius = 0.3;
    pressure = (*q_inv_pt) * 2.0 / (*alpha_pt) *
               (1.0 + 1.0 / (*alpha_pt) / bubble_radius);
  }
} // namespace relaxing_bubble

#endif
