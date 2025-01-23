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
    double wall_velocity = 1.0;
    int max_adapt = 0;
    std::string output_directory = "RESLT";
    std::string restart_filename = "";
    unsigned azimuthal_mode_number = 0;
    unsigned bulk_element_number_of_plot_points = 3;
    unsigned error_estimator_flag = 1;
    unsigned initial_number_of_free_surface_points = 32;
    unsigned interval_between_adapts = 0;
    unsigned max_newton_iterations = 40;
    unsigned max_number_of_adapts_for_refinement = 0;
    unsigned surface_element_number_of_plot_points = 3;
  };

  Params create_parameters_from_file(const std::string& filename)
  {
    Params params;

    std::ifstream parameter_filestream(filename);
    std::string input_string;

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    *params.reynolds_inverse_froude_number_pt = stod(input_string);

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
    input_string.erase(
      remove_if(input_string.begin(), input_string.end(), isspace),
      input_string.end());
    if (input_string.size() > 0)
    {
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

    parameter_filestream.close();

    return params;
  }

  void save_parameters_to_file(const Params& params,
                               const std::string& filename)
  {
    std::ofstream parameter_filestream(filename);

    parameter_filestream << std::setprecision(
      std::numeric_limits<double>::max_digits10);

    parameter_filestream << *params.reynolds_inverse_froude_number_pt
                         << " # Bond number" << "\n";
    parameter_filestream << params.capillary_number << " # Capillary number"
                         << "\n";
    parameter_filestream << params.reynolds_number << " # Reynolds number"
                         << "\n";
    parameter_filestream << params.contact_angle * 180.0 /
                              oomph::MathematicalConstants::Pi
                         << " # Contact angle" << "\n";
    parameter_filestream << params.max_adapt << " # Max number of adapt steps"
                         << "\n";
    parameter_filestream << params.output_directory << " # Output directory"
                         << "\n";
    parameter_filestream << params.slip_length << " # Slip length" << "\n";
    parameter_filestream << params.min_element_length
                         << " # Mininum element length" << "\n";
    parameter_filestream << params.element_length_ratio
                         << " # Element length ratio" << "\n";
    parameter_filestream << params.flux_withdraw_speed << " # Withdraw speed "
                         << "\n";
    parameter_filestream << params.flux_duration << " # Flux duration" << "\n";
    parameter_filestream << params.ramp_up_time << " # Flux ramp up time"
                         << "\n";
    parameter_filestream << params.max_element_size << " # Max element area"
                         << "\n";
    parameter_filestream << params.min_element_size << " # Min element area"
                         << "\n";
    parameter_filestream << params.max_free_surface_polyline_length
                         << " # Max_free_surface_polyline_length " << "\n";
    parameter_filestream << params.max_slip_polyline_length
                         << " # Max_slip_polyline_length" << "\n";
    parameter_filestream << params.interval_between_adapts
                         << " # interval_between_adapts" << "\n";
    parameter_filestream
      << params.error_estimator_flag
      << " # Error estimator flag, 0 ContactLine, 1 Z2, 2 Corner" << "\n";
    parameter_filestream << params.nu << " # Pseudo-solid Poisson ratio (Nu)"
                         << "\n";
    parameter_filestream << params.temporal_tolerance << " # Temporal tolerance"
                         << "\n";
    parameter_filestream << params.max_timestep << " # Max timestep" << "\n";
    parameter_filestream << params.is_adaptive_timestepping
                         << " # Use adaptive timestepping" << "\n";
    parameter_filestream << params.max_permitted_z2_error
                         << " # Max permitted Z2 error" << "\n";
    parameter_filestream << params.min_permitted_z2_error
                         << " # Min permitted Z2 error" << "\n";
    parameter_filestream << params.final_time << " # Target final time" << "\n";
    parameter_filestream << params.time_step << " # Time step" << "\n";
    parameter_filestream << params.restart_filename << " # Restart filename"
                         << "\n";
    parameter_filestream << params.wall_velocity << " # Wall velocity" << "\n";
    parameter_filestream << params.azimuthal_mode_number
                         << " # Azimuthal mode number" << "\n";
    parameter_filestream
      << params.max_number_of_adapts_for_refinement
      << " # Max number of adapts for initial mesh refinement" << "\n";
    parameter_filestream << params.is_strong_contact_angle
                         << " # Use strong contact angle" << "\n";
    parameter_filestream << params.polyline_refinement_tolerence
                         << " # Free surface polyline refinement tolerence"
                         << "\n ";
    parameter_filestream << params.augmented_radius
                         << " # Augmented region's radius" << std::endl;

    parameter_filestream.close();
  }
} // namespace oomph
#endif
