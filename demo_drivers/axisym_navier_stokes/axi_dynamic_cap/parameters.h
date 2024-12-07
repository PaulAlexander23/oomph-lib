#ifndef PARAMETERS_HEADER
#define PARAMETERS_HEADER

#include <sys/stat.h>
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
    double newton_solver_tolerance = 1e-8;
    double nu = 0.25;
    double polyline_refinement_tolerence = 4e-3; // 8e-3
    double polyline_unrefinement_tolerence = 2e-3; // 4e-3
    double ramp_up_time = 0.1;
    double reynolds_inverse_froude_number = 0.0;
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
    double wall_velocity = 1.0;
    int max_adapt = 0;
    std::string output_directory = "";
    std::string restart_filename = "";
    unsigned azimuthal_mode_number = 0;
    unsigned bulk_element_number_of_plot_points = 3;
    unsigned error_estimator_flag = 1;
    unsigned initial_number_of_free_surface_points = 32;
    unsigned interval_between_adapts = 5;
    unsigned max_newton_iterations = 20;
    unsigned max_number_of_adapts_for_refinement = 20;
    unsigned surface_element_number_of_plot_points = 3;
  };

  Params create_parameters_from_file(const std::string& filename)
  {
    Params params;

    std::ifstream parameter_filestream(filename);
    std::string input_string;

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.reynolds_inverse_froude_number = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.capillary_number = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.reynolds_number = stod(input_string);
    params.reynolds_strouhal_number =
      params.reynolds_number * params.strouhal_number;

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.contact_angle =
      stod(input_string) * oomph::MathematicalConstants::Pi / 180.0;

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.max_adapt = stoi(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    if (input_string.size() > 0) input_string.resize(input_string.size() - 1);
    params.output_directory = input_string;

    struct stat stat_buffer;
    if (stat(params.output_directory.c_str(), &stat_buffer))
    {
      oomph_info << "WARNING: Directory doesn't exist." << std::endl;
    }

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.slip_length = std::stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.min_element_length = stod(input_string);
    params.inner_min_element_length = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.element_length_ratio = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.flux_withdraw_speed = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.flux_duration = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.ramp_up_time = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.max_element_size = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.min_element_size = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.max_free_surface_polyline_length = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.max_slip_polyline_length = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.interval_between_adapts = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.error_estimator_flag = stoi(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.nu = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.temporal_tolerance = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.max_timestep = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.is_adaptive_timestepping = stoi(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.max_permitted_z2_error = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.min_permitted_z2_error = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.final_time = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.time_step = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    if (input_string.size() > 0)
    {
      input_string.resize(input_string.size() - 1);
      params.is_restarting = true;
      params.restart_filename = input_string;
    }

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.wall_velocity = stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.azimuthal_mode_number = stoi(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.max_number_of_adapts_for_refinement = stoi(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.is_strong_contact_angle = stoi(input_string);

    try
    {
      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      params.polyline_refinement_tolerence = stod(input_string);

      params.polyline_unrefinement_tolerence =
        0.5 * params.polyline_refinement_tolerence;
    }
    catch (std::exception& e)
    {
    }

    try
    {
      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      params.augmented_radius = stod(input_string);
    }
    catch (std::exception& e)
    {
    }

    return params;
  }
} // namespace oomph
#endif
