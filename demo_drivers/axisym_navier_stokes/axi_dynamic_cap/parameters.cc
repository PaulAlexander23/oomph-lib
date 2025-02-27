#include "parameters.h"

#include <ctype.h>
#include <sys/stat.h>
#include <algorithm>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>

#include "generic/oomph_definitions.h"

namespace oomph{
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
    *params.wall_velocity_pt = stod(input_string);

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

    try
    {
      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      params.initial_number_of_free_surface_points = stoi(input_string);
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
    parameter_filestream << *params.wall_velocity_pt << " # Wall velocity"
                         << "\n";
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
                         << " # Augmented region's radius" << "\n";
    parameter_filestream << params.initial_number_of_free_surface_points
                         << " # Initial number of free surface points";
    parameter_filestream << std::endl;
    parameter_filestream.close();
  }
};
