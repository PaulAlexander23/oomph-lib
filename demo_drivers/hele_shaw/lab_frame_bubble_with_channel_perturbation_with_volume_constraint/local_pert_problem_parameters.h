#ifndef PROBLEM_PARAMETER_HEADER
#define PROBLEM_PARAMETER_HEADER

#include "generic.h"
//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
namespace Problem_Parameter
{
  /// Doc info object
  DocInfo Doc_info;

  /// Restart input filename
  std::string restart_input_filename;

  /// Pseudo-solid Poisson ratio
  double Nu = 0.3;

  // double Radius=0.7;


  bool Minisolve_state = false;
  bool ignore_height_effects_in_dynamic_bc = false;
  bool ignore_height_effects_in_viscous_terms = false;

  bool Set_steady_conditions_after_adapt = false;

  double* global_G_pt = 0;
  double* global_G_upper_outflow_pt = 0;
  double* global_Q_inv_pt = 0;
  double* global_Frame_speed_pt = 0;
  double* global_Obstacle_height_pt = 0;
  double* global_Obstacle_width_pt = 0;
  double* global_alpha_pt = 0;
  double* global_Asymmetry_pt = 0;
  double* global_Bubble_pressure_pt = 0;
  double* global_St_pt = 0;
  double* global_delta_pt = 0;
  double* global_CoM_X_pt = 0;
  double* global_CoM_Y_pt = 0;
  double* global_measure_pt = 0;


  double Major_Radius;
  double Minor_Radius;
  double xcenter;
  double ycenter;
  unsigned circpts;

  unsigned n_integral_measures = 12;

  // Node* Tip_node_pt;
  bool ignore_bulk_effects = false;

  /// \short Volume of the bubble (negative because it's outside the
  /// fluid!)
  double Volume = -MathematicalConstants::Pi * Major_Radius * Minor_Radius;
  double Centre_of_mass = 0;

  /// \short Scaling factor for inflow velocity (allows it to be switched off
  /// to do hydrostatics)
  double Inflow_veloc_magnitude = 0.0;

  /// \short Length of the channel
  double Length = 4.0;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw* Constitutive_law_pt = 0;

  /// Trace file
  ofstream Trace_file;

  /// \short File to document the norm of the solution (for validation
  /// purposes -- triangle doesn't give fully reproducible results so
  /// mesh generation/adaptation may generate slightly different numbers
  /// of elements on different machines!)
  ofstream Norm_file;

  ofstream OccluHeight_file;

  void channel_height_function(const Vector<double>& x, double& b)
  {
    /// This function should have obstacle width and height, asymmetry and
    /// possibly sharpness.

    double height = 0.0;
    double rms_width = 0.1;
    double centre_x = 0.0;
    double centre_y = 0.0;

    if (global_Obstacle_height_pt != 0)
    {
      height = *global_Obstacle_height_pt;
    }
    if (global_Obstacle_width_pt != 0)
    {
      rms_width = *global_Obstacle_width_pt;
    }

    // Transform the domain as required i.e. from -1:1 to 0:1
    double local_x = x[0];
    double local_y = x[1];

    b = 1.0 - height * std::exp(-(local_x - centre_x) * (local_x - centre_x) /
                                  (2 * rms_width * rms_width) -
                                (local_y - centre_y) * (local_y - centre_y) /
                                  (2 * rms_width * rms_width));
  }

  ofstream UpperWall_file;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    channel_height_function(x, b);
    if (ignore_bulk_effects)
    {
      b = 1;
      UpperWall_file << b << std::endl;
    }
    dbdt = 0;
  }

  void frame_speed_fct(const Vector<double>& x, Vector<double>& U_frame)
  {
    if (global_Frame_speed_pt != 0)
    {
      U_frame[0] = *global_Frame_speed_pt;
      U_frame[1] = 0;
    }
    else
    {
      U_frame[0] = 0;
      U_frame[1] = 0;
    }
  }

  void bubble_pressure_function(const Vector<double>& x, double& pressure)
  {
    /// For slight simplicity, we can put the transverse curvature term as a
    /// spatially varying bubble pressure.

    double height = 1.0;
    if (ignore_height_effects_in_dynamic_bc == false)
    {
      channel_height_function(x, height);
    }
    double alpha = *global_alpha_pt;
    double Q = 1 / (*global_Q_inv_pt);
    pressure = -1.0 / (3.0 * Q * alpha * height) + (*global_Bubble_pressure_pt);
  }

  void bubble_pressure_gradient_function(const Vector<double>& x,
                                         Vector<double>& grad_pb)
  {
    /// ALICE_FLAG
    std::cout << "Evaluating bubble pressure gradient" << std::endl;
    std::cout << "Missing height dependence!" << std::endl;
  }

  void normal_flux_ahead_of_bubble(const Vector<double>& x, double& flux)
  {
    /// Ahead of the bubble we have the boundary condition dp/dx = -G.
    /// We need to supply the function b^3 n.grad p = -G b^3.
    double G;
    if (global_G_pt != 0)
    {
      G = *global_G_pt;
    }
    else
    {
      G = 1;
    }
    double p_x = -1 * G;
    flux = p_x;
  }

  void normal_flux_behind_bubble(const Vector<double>& x, double& flux)
  {
    /// Ahead of the bubble we have the boundary condition dp/dx = -G.
    /// We need to supply the function b^3 n.grad p = -G b^3.
    double G;
    if (global_G_pt != 0)
    {
      G = *global_G_pt;
    }
    else
    {
      G = 1;
    }
    double p_x = G;
    flux = p_x;
  }

  Vector<double> vector_of_eigenvalues_rp(10, 2.0);
  Vector<double> vector_of_eigenvalues_ip(10, 2.0);

  ofstream M_file;
} // namespace Problem_Parameter
#endif
