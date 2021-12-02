#ifndef OOMPH_HELE_SHAW_PARAMETER_CONVERTER
#define OOMPH_HELE_SHAW_PARAMETER_CONVERTER

namespace oomph
{
  void convert_to_capillary_nondimensionalisation(const double& r,
                            const double& width,
                            const double& volume,
                            const double& Q,
                            double& length_ratio,
                            double& pressure_ratio,
                            double& velocity_ratio,
                            double& time_ratio,
                            double& new_r,
                            double& new_width,
                            double& new_volume,
                            double& new_Ca)
  {
    /// Ratios of old to new
    length_ratio = 0.5;
    pressure_ratio = 12;
    velocity_ratio = length_ratio;
    time_ratio = 1;

    /// New parameters
    new_r = r * length_ratio;
    new_width = width * length_ratio;
    new_volume = volume * length_ratio * length_ratio;
    new_Ca = Q * velocity_ratio;
  }

} // namespace oomph

#endif
