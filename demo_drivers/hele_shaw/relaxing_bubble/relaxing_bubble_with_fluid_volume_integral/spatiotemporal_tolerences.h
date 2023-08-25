#ifndef OOMPH_SPATIOTEMPORAL_TOLERANCES_HEADER
#define OOMPH_SPATIOTEMPORAL_TOLERANCES_HEADER

#include <iostream>
#include "generic.h"

namespace oomph
{
  class SpatiotemporalTolerances
  {
  public:
    SpatiotemporalTolerances()
    {
      // --- Default values ---
      // Polyline arguments
      polyline_refinement_tolerance = 0.08;
      polyline_unrefinement_tolerance = 0.04;
      maximum_polyline_segment_length = 0.02;
      initial_number_of_polynomial_vertices = 64;

      // Triangle bulk mesh arguments
      maximum_element_size = 5e-2; // 5e-2
      minimum_element_size = 1e-6; // 1e-6
      maximum_permitted_error = 2e-5; // 5e-6; // 2e-5
      minimum_permitted_error = 5e-6; //1e-6; // 5e-6
      initial_target_element_area = 1e-1;

      // Timestepper??
      adaptive_timestepping = true;
      initial_timestep = 1e-2;

      // Remeshing
      remesh_interval = 5;
      remesh_initial_condition = false;
    }

    void doc(DocInfo& doc_info, string filename)
    {
      string full_filename = doc_info.directory() + filename;

      ofstream output_stream;
      output_stream.open(full_filename);

      output_stream << "polyline_refinement_tolerance: "
                    << polyline_refinement_tolerance << "\n";
      output_stream << "polyline_unrefinement_tolerance: "
                    << polyline_unrefinement_tolerance << "\n";
      output_stream << "maximum_polyline_segment_length: "
                    << maximum_polyline_segment_length << "\n";
      output_stream << "initial_number_of_polynomial_vertices: "
                    << initial_number_of_polynomial_vertices << "\n";
      output_stream << "maximum_element_size: " << maximum_element_size << "\n";
      output_stream << "minimum_element_size: " << minimum_element_size << "\n";
      output_stream << "maximum_permitted_error: " << maximum_permitted_error
                    << "\n";
      output_stream << "minimum_permitted_error: " << minimum_permitted_error
                    << "\n";
      output_stream << "initial_target_element_area: "
                    << initial_target_element_area << "\n";
      output_stream << "adaptive_timestepping: " << adaptive_timestepping
                    << "\n";
      output_stream << "initial_timestep: " << initial_timestep << "\n";
      output_stream << "remesh_interval: " << remesh_interval << "\n";
      output_stream << "remesh_initial_condition: " << remesh_initial_condition
                    << "\n";

      output_stream.close();
    }

    void set_polyline_refinement_tolerance(double new_value)
    {
      polyline_refinement_tolerance = new_value;
    }

    double get_polyline_refinement_tolerance()
    {
      return polyline_refinement_tolerance;
    }

    void set_polyline_unrefinement_tolerance(double new_value)
    {
      polyline_unrefinement_tolerance = new_value;
    }

    double get_polyline_unrefinement_tolerance()
    {
      return polyline_unrefinement_tolerance;
    }

    void set_maximum_polyline_segment_length(double new_value)
    {
      maximum_polyline_segment_length = new_value;
    }

    double get_maximum_polyline_segment_length()
    {
      return maximum_polyline_segment_length;
    }

    void set_initial_number_of_polynomial_vertices(unsigned new_value)
    {
      initial_number_of_polynomial_vertices = new_value;
    }

    unsigned get_initial_number_of_polynomial_vertices()
    {
      return initial_number_of_polynomial_vertices;
    }

    void set_maximum_element_size(double new_value)
    {
      maximum_element_size = new_value;
    }

    double get_maximum_element_size()
    {
      return maximum_element_size;
    }

    void set_minimum_element_size(double new_value)
    {
      minimum_element_size = new_value;
    }

    double get_minimum_element_size()
    {
      return minimum_element_size;
    }

    void set_maximum_permitted_error(double new_value)
    {
      maximum_permitted_error = new_value;
    }

    double get_maximum_permitted_error()
    {
      return maximum_permitted_error;
    }

    void set_minimum_permitted_error(double new_value)
    {
      minimum_permitted_error = new_value;
    }

    double get_minimum_permitted_error()
    {
      return minimum_permitted_error;
    }

    void set_initial_target_element_area(double new_value)
    {
      initial_target_element_area = new_value;
    }

    double get_initial_target_element_area()
    {
      return initial_target_element_area;
    }

    void set_adaptive_timestepping(bool new_value)
    {
      adaptive_timestepping = new_value;
    }

    bool get_adaptive_timestepping()
    {
      return adaptive_timestepping;
    }

    void set_initial_timestep(double new_value)
    {
      initial_timestep = new_value;
    }

    double get_initial_timestep()
    {
      return initial_timestep;
    }

    void set_remesh_interval(unsigned new_value)
    {
      remesh_interval = new_value;
    }

    unsigned get_remesh_interval()
    {
      return remesh_interval;
    }

    void set_remesh_initial_condition(bool new_value)
    {
      remesh_initial_condition = new_value;
    }

    bool get_remesh_initial_condition()
    {
      return remesh_initial_condition;
    }

  private:
    // Polyline arguments
    double polyline_refinement_tolerance;
    double polyline_unrefinement_tolerance;
    double maximum_polyline_segment_length;
    unsigned initial_number_of_polynomial_vertices;

    // Triangle bulk mesh arguments
    double maximum_element_size;
    double minimum_element_size;
    double maximum_permitted_error;
    double minimum_permitted_error;
    double initial_target_element_area;

    // Timestepper??
    bool adaptive_timestepping;
    double initial_timestep;

    // Remeshing
    unsigned remesh_interval;
    bool remesh_initial_condition;
  };


}; // namespace oomph
#endif
