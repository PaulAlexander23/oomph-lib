#include <iostream>
#include <chrono>

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

#include "relaxing_bubble_problem.h"
#include "integral_problem.h"

using namespace oomph;
using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])
{
  /// Start timer
  auto start_time = high_resolution_clock::now();

  /// State program name
  cout << "Moving bubble with total volume integral" << endl;

  /// Default values for the command line arguments
  double Q = 0.05;
  string output_directory = "RESLT/";
  bool has_unrecognised_arg = false;

  /// Setup and parse command line arguments
  CommandLineArgs::setup(argc, argv);
  /// OPTIONAL
  CommandLineArgs::specify_command_line_flag(
    "-q", &Q, "Optional: Capillary/Channel flux parameter. Default = 0.05");
  CommandLineArgs::specify_command_line_flag(
    "-o", &output_directory, "Optional: Output directory. Default = RESLT/");
  CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);

  /// Typically varied parameters
  const double circular_radius = 0.3;
  const double volume =
    MathematicalConstants::Pi * circular_radius * circular_radius;
  const double h = 0.0;

  /// Set other fixed parameters
  const double major_radius = 0.3;
  const double width = 0.25;

  /// Convert into the Capillary scaling
  double Ca;
  double new_volume;
  double new_h;
  convert_parameters_from_q_nd_to_capillary_nd(
    Q, volume, h, Ca, new_volume, new_h);

  double length_ratio;
  double depth_ratio;
  double time_ratio;
  double pressure_ratio;
  get_ratios_of_q_nd_to_capillary_nd(
    length_ratio, depth_ratio, time_ratio, pressure_ratio);

  relaxing_bubble::total_flux = 0.0;
  relaxing_bubble::alpha = 40.0;
  relaxing_bubble::major_radius = major_radius * length_ratio;
  relaxing_bubble::ca_inv = 1 / Ca;
  relaxing_bubble::st = 1.0;
  relaxing_bubble::bubble_initial_centre_y = length_ratio * (1 + 0.01);
  relaxing_bubble::perturbation_amplitude = new_h;
  relaxing_bubble::perturbation_rms_width = length_ratio * width;
  relaxing_bubble::target_bubble_volume = new_volume;
  relaxing_bubble::nu = 0.3;
  // Create generalised Hookean constitutive equations
  relaxing_bubble::constitutive_law_pt =
    new GeneralisedHookean(&relaxing_bubble::nu);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;
  doc_info.set_directory(output_directory);

  /// Compute total volume of channel
  IntegralProblem<QIntegralElement<3>> integral_problem(
    relaxing_bubble::channel_depth);
  integral_problem.newton_solve();
  integral_problem.doc_solution(doc_info);

  relaxing_bubble::total_volume = integral_problem.result();
  cout << "Total volume: " << relaxing_bubble::total_volume << endl;

  /// Set target fluid volume
  relaxing_bubble::target_fluid_volume =
    (relaxing_bubble::total_volume) - (relaxing_bubble::target_bubble_volume);

  /// Print the parameters
  relaxing_bubble::print_parameters();
  relaxing_bubble::doc_parameters(doc_info, "parameters.dat");

  SpatiotemporalTolerances* tolerances_pt = new SpatiotemporalTolerances;

  /// I think the polyline segment length is in terms of the local coordinate
  /// which doesn't require the scaling transform.
  // tolerances_pt->set_maximum_polyline_segment_length(
  //  tolerances_pt->get_maximum_polyline_segment_length() * length_ratio);

  tolerances_pt->set_maximum_element_size(
    tolerances_pt->get_maximum_element_size() * length_ratio * length_ratio);
  tolerances_pt->set_minimum_element_size(
    tolerances_pt->get_minimum_element_size() * length_ratio * length_ratio);
  tolerances_pt->set_initial_target_element_area(
    tolerances_pt->get_initial_target_element_area() * length_ratio *
    length_ratio);

  tolerances_pt->doc(doc_info, "tolerances.dat");

  /// Create problem
  RelaxingBubbleProblem<MyNewElementWithIntegral> problem(tolerances_pt);

  /// Run self tests
  bool run_self_test = false;
  if (run_self_test)
  {
    if (problem.self_test())
    {
      throw OomphLibError(
        "Self test failed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  /// Solve for initial conditions
  problem.solve_for_initial_conditions(doc_info);

  /// Set initial timestep and simulation duration
  double dt = 1e-2 * time_ratio;
  double tF = 4e0 * time_ratio;

  /// Iterate the adaptive step timestepper until the final time, documenting
  /// the solution
  problem.iterate_timestepper(dt, tF, doc_info);

  /// Delete the constitutive law
  delete relaxing_bubble::constitutive_law_pt;

  /// End timer
  auto end_time = high_resolution_clock::now();
  duration<double> diff = duration_cast<minutes>(end_time - start_time);
  cout << "Total real-time duration: " << diff.count() << endl;

  return 0;
}
