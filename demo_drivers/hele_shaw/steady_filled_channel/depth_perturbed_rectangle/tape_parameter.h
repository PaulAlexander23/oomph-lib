#ifndef OOMPH_STATIC_HS_PROBLEM_PARAMETER_HEADER
#define OOMPH_STATIC_HS_PROBLEM_PARAMETER_HEADER

#include "generic.h"

namespace problem_parameter
{
  const double domain_length = 1.0;
  const double domain_width = 1.0;

  double* inlet_area_pt = 0;
  double* outlet_area_pt = 0;

  /// This is non-dimensionalised to 1
  const double total_flux = 1.0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    double tape_height = 0.1;
    double tape_width = 0.4;
    double tape_sharpness = 40;
    double tape_centre_y = 0.5;

    b = 1.0 - tape_height * 0.5 *
                (std::tanh(tape_sharpness * (x[1] - tape_centre_y + 0.5 * tape_width)) -
                 std::tanh(tape_sharpness * (x[1] - tape_centre_y - 0.5 * tape_width)));
    dbdt = 0.0;
  }

  void get_dirichlet_bc(const Vector<double>& x, double& p)
  {
    /// At the outlet we set the pressure to be zero
    p = 0.0;
  }

  void get_inlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_area and total flux

    double dpdx = total_flux / *inlet_area_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * pow(b, 3.0) / 12.0;
  }

  void get_outlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_area and total flux

    double dpdx = total_flux / *outlet_area_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * pow(b, 3.0) / 12.0;
  }

} // namespace problem_parameter
#endif
