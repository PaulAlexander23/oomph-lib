#ifndef OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER
#define OOMPH_RELAXING_BUBBLE_PARAMETERS_HEADER

#include "generic.h"

namespace moving_bubble
{
  /// Dimensionless parameters pointer
  double* ca_inv_pt = 0;
  double* st_pt = 0;
  double* alpha_pt = 0;
  const double total_flux = 1.0;

  double major_radius;

  /// Pseudo-solid Poisson ratio pointer
  double* nu_pt = 0;
  /// Constitutive law pointer
  ConstitutiveLaw* constitutive_law_pt = 0;

  /// Bubble pressure pointer
  double* bubble_pressure_pt = 0;

  /// Volume pointers
  double* target_bubble_volume_pt = 0;
  double* target_fluid_volume_pt = 0;
  double* total_volume_pt = 0;

  /// Flux area pointers
  double* inlet_area_pt = 0;
  double* outlet_area_pt = 0;

  /// Channel depth without rate of change
  void channel_depth(const Vector<double>& x, double& b)
  {
    double amplitude = 0.024;
    double rms_width = 0.1;
    double centre_x = 1.5;
    double centre_y = 0.5;

    // Transform the domain as required i.e. from -1:1 to 0:1
    double local_x = x[0];
    double local_y = x[1];

    b = 1.0 - amplitude * exp(-(local_x - centre_x) * (local_x - centre_x) /
                                (2 * rms_width * rms_width) -
                              (local_y - centre_y) * (local_y - centre_y) /
                                (2 * rms_width * rms_width));
  }

  /// Channel depth with rate of change
  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    channel_depth(x, b);
    dbdt = 0.0;
  }

  // Wall speed function or moving reference frame speed function
  void wall_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    U_frame[0] = 0.0;
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
    /// upper wall function, inlet_area and total flux

    double dpdx = total_flux / *inlet_area_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b);
  }

  /// Outlet boundary condition
  void get_outlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_area and total flux
    double dpdx = -total_flux / *outlet_area_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b);
  }

} // namespace moving_bubble

#endif
