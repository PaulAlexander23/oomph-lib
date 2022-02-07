#ifndef OOMPH_HELE_SHAW_PARAMETER_CONVERTER
#define OOMPH_HELE_SHAW_PARAMETER_CONVERTER

namespace oomph
{
  /// Get the ratios of the different independent variables scaled first by
  /// fixing the pressure gradient at the end of the domain to 1, to the
  /// second method of fixing the flux at the end of the channel to 1.
  void get_ratios_of_q_nd_to_capillary_nd(double& length_ratio,
                                          double& depth_ratio,
                                          double& time_ratio,
                                          double& pressure_ratio)
  {
    /// I = int_{-1}^{1} (\bar{b}(y))^3 dy
    /// where b =
    /// Calculated using a simple finite element code.
    /// Using 2 elements
    // const double I = 1.9648098177825899;
    /// Using 40 elements
    // const double I = 1.9648144969660464;
    /// Using 80 elements, accurate to at least 6dp
    const double I = 1.9648144083521897;

    /// Ratios of old to new
    length_ratio = 0.5;
    depth_ratio = 1.0;
    time_ratio = I / 4.0;
    pressure_ratio = 1.0 / I; // or 12.0 / I;??
  }

  /// Convert the parameters used from a nondimensionalation which
  /// fixes the pressure gradient at the end of the domain to 1, to a
  /// second method which fixes the flux at the end of the channel to 1.
  void convert_parameters_from_q_nd_to_capillary_nd(const double& Q,
                                                    const double& V,
                                                    const double& h,
                                                    double& Ca,
                                                    double& new_V,
                                                    double& new_h)
  {
    double length_ratio;
    double depth_ratio;
    double time_ratio;
    double pressure_ratio;
    get_ratios_of_q_nd_to_capillary_nd(
      length_ratio, depth_ratio, time_ratio, pressure_ratio);
    Ca = 6.0 * Q / pressure_ratio;
    new_V = length_ratio * length_ratio * depth_ratio * V;
    new_h = depth_ratio * h;
  }

  /// Convert the parameters used from a nondimensionalation which
  /// fixes the pressure gradient at the end of the domain to 1, to a
  /// second method which fixes the flux at the end of the channel to 1.
  void convert_parameters_from_capillary_nd_to_q_nd(const double& Ca,
                                                    const double& V,
                                                    const double& h,
                                                    double& Q,
                                                    double& new_V,
                                                    double& new_h)
  {
    double length_ratio;
    double depth_ratio;
    double time_ratio;
    double pressure_ratio;
    get_ratios_of_q_nd_to_capillary_nd(
      length_ratio, depth_ratio, time_ratio, pressure_ratio);
    Q = pressure_ratio / 6.0 * Ca;
    new_V = V / (length_ratio * length_ratio * depth_ratio);
    new_h = h / depth_ratio;
  }

} // namespace oomph

#endif
