#ifndef OOMPH_STATIC_HS_PROBLEM_PARAMETER_HEADER
#define OOMPH_STATIC_HS_PROBLEM_PARAMETER_HEADER

#include "generic.h"

namespace problem_parameter
{
  const double domain_length = 2.0;
  const double domain_width = 1.0;

  double* inlet_integral_pt = 0;
  double* outlet_integral_pt = 0;

  /// This is non-dimensionalised to 1
  const double total_flux = 1.0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    b = 1.0 - 0.2 * std::exp(-pow(x[0] - 1, 2.0) / (2.0 * pow(0.1, 2.0)) -
                             pow(x[1] - 0.5, 2.0) / (2.0 * pow(0.1, 2.0)));
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
    /// upper wall function, inlet_integral and total flux

    double dpdx = total_flux / *inlet_integral_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * pow(b, 3.0) / 12.0;
  }

  void get_outlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_integral and total flux

    double dpdx = total_flux / *outlet_integral_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * pow(b, 3.0) / 12.0;
  }

} // namespace problem_parameter
#endif
