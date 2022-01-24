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
  double* bubble_pressure_pt = 0;
  double* target_bubble_volume_pt = 0;
  double* target_fluid_volume_pt = 0;
  double* total_volume_pt = 0;

  ConstitutiveLaw* constitutive_law_pt = 0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    double amplitude = 0.2;
    double sharpness = 40;
    double width = 0.25;
    double x_centre = 0.5;
    double y_centre = 0.125;

    // b = 1.0 - amplitude / 2 *
    //            ((tanh(sharpness * (x[0] - x_centre + width / 2)) -
    //              tanh(sharpness * (x[0] - x_centre - width / 2))) /
    //               2 +
    //             (tanh(sharpness * (x[1] - y_centre + width / 2)) -
    //              tanh(sharpness * (x[1] - y_centre - width / 2))) /
    //               2);

    b =
      1.0 - amplitude / 2 *
              (cos(2 * MathematicalConstants::Pi * (x[0] - x_centre) / width) +
               cos(2 * MathematicalConstants::Pi * (x[1] - y_centre) / width));

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
