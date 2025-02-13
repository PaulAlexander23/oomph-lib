#ifndef PARAMETERS_HEADER
#define PARAMETERS_HEADER

#include <sys/stat.h>
#include <limits>
#include <memory>

#include "generic.h"

namespace oomph
{
  struct Params
  {
    Vector<double> gravity_vector = {0, -1.0, 0};
    bool is_adaptive_timestepping = false;
    bool is_restarting = false;
    bool is_strong_contact_angle = false;
    double augmented_radius = 0.3;
    double capillary_number = 1.0;
    double contact_angle = 2.0 / 3.0 * MathematicalConstants::Pi;
    double element_length_ratio = 1.5;
    double final_time = 0.0;
    double flux_duration = 1e0;
    double flux_withdraw_speed = 0e0;
    double free_surface_error_tolerence = 8e-3;
    double initial_fluid_height = 3.5;
    double inner_min_element_length = 2e-4;
    double max_element_size = 0.5 * std::pow(5e-1, 2.0);
    double max_free_surface_polyline_length = 5e-2;
    double max_permitted_z2_error = 1e-3;
    double max_residual = 1e3;
    double max_slip_polyline_length = 1e-1;
    double max_timestep = 1e-1;
    double min_element_length = 2e-4;
    double min_element_size = 0.5 * std::pow(min_element_length, 2.0);
    double min_permitted_angle = 15;
    double min_permitted_z2_error = 1e-7;
    double min_permitted_mesh_residual = 0.01;
    double max_permitted_mesh_residual = 1.0;
    double newton_solver_tolerance = 1e-8;
    double nu = 0.25;
    double polyline_refinement_tolerence = 4e-3; // 8e-3
    double polyline_unrefinement_tolerence = 2e-3; // 4e-3
    double ramp_up_time = 0.1;
    double* reynolds_inverse_froude_number_pt = new double(0.0);
    double reynolds_number = 0.0;
    double reynolds_strouhal_number = 0.0;
    double right_angle = MathematicalConstants::Pi * 90.0 / 180.0;
    double sigma = 0.0;
    double slip_length = 1e0;
    double small_r = 1e-4;
    double strouhal_number = 1.0;
    double temporal_tolerance = 1e0;
    double time_step = 0.0;
    double uniform_element_area = 0.5 * std::pow(5e-1, 2.0);
    double viscosity_ratio = 1.0;
    double volume = 3.5 / 2.0;
    double* wall_velocity_pt = new double(0.0);
    int max_adapt = 0;
    std::string output_directory = "RESLT";
    std::string restart_filename = "";
    int azimuthal_mode_number = 0;
    unsigned bulk_element_number_of_plot_points = 3;
    unsigned error_estimator_flag = 1;
    unsigned initial_number_of_free_surface_points = 32;
    unsigned interval_between_adapts = 0;
    unsigned max_newton_iterations = 40;
    unsigned max_number_of_adapts_for_refinement = 0;
    unsigned surface_element_number_of_plot_points = 3;
  };

  Params create_parameters_from_file(const std::string& filename);

  void save_parameters_to_file(const Params& params,
                               const std::string& filename);
} // namespace oomph
#endif
