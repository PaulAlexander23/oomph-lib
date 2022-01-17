#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "ode.h"

using namespace oomph;
using namespace std;

namespace problem_parameter
{
  /// Global data
  double* global_flux_pt = 0;
  double* global_frame_speed_pt = 0;
  double* global_frame_travel_pt = 0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    /// This function should have obstacle width and height, asymmetry and
    /// possibly sharpness.

    double height = 0.1;
    double rms_width = 0.1;
    double centre_x = 1 + *global_frame_travel_pt;
    double centre_y = 0.5;

    // Transform y such that the domain is between 0 and 1 rather than -1 and 1
    double local_x = x[0];
    double local_y = x[1];
    double exp_term =
      height * exp(-(local_x - centre_x) * (local_x - centre_x) /
                     (2 * rms_width * rms_width) -
                   (local_y - centre_y) * (local_y - centre_y) /
                     (2 * rms_width * rms_width));

    b = 1.0 - exp_term;

    dbdt = *global_frame_speed_pt * 2.0 * (local_x - centre_x) * exp_term;
  }

  void get_dirichlet_bc(const Vector<double>& x, double& p)
  {
    /// At the outlet we set the pressure to be zero
    p = 0.0;
  }

  void get_neumann_bc(const Vector<double>& x, double& dpdx)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_area and total flux
    double b = 0.0;
    double dbdt = 0.0;
    upper_wall_fct(x, b, dbdt);
    double total_flux = 1.0;
    double inlet_area = 1.0;
    dpdx = total_flux / inlet_area * (b * b * b);
  }

  class ODEFunctor : public SolutionFunctorBase
  {
  public:
    /// Constructor
    ODEFunctor() {}

    /// Destructor
    virtual ~ODEFunctor() {}

    /// Exact or approximate solution. Used for initialisation and/or error
    /// checking
    Vector<double> operator()(const double& t, const Vector<double>& x) const
    {
      Vector<double> output(1);

      output[0] = sin(4 * atan(1) * t);

      return output;
    }

    /// Derivative function. Specifies the ODE that we are solving
    Vector<double> derivative(const double& t,
                              const Vector<double>& x,
                              const Vector<double>& u) const
    {
      Vector<double> output(1);

      output[0] = 4 * atan(1) * cos(4 * atan(1) * t);

      return output;
    }
  };
} // namespace problem_parameter

template<class ELEMENT>
class HeleShawChannelProblem : public Problem
{
public:
  /// Constructor
  HeleShawChannelProblem();

  /// Destructor (empty)
  ~HeleShawChannelProblem() {}

  /// Set the initial conditions including previous history
  void set_initial_condition(const double& t_step, DocInfo& doc_info);

  /// Iterate forward in time in steps of t_step until t_final
  void iterate_timestepper(const double& t_step,
                           const double& t_final,
                           DocInfo& doc_info);

  /// Doc the solution
  void doc_solution(DocInfo& doc_info);

private:
  /// Generate data
  void generate_data();

  /// Generate mesh
  void generate_mesh();

  /// Create flux elements
  void create_flux_elements(const unsigned& boundary);

  /// Generate mesh
  void assign_mesh();

  /// Pin dirichlet outlet boundary
  void pin_dirichlet_boundaries();

  /// Upcast elements and finalise setup
  void setup_elements();

  /// Set boundary condition values
  void set_boundary_conditions();

  /// Update the problem specs before solve
  void actions_before_newton_solve();

  /// Update the problem specs before solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before timestep
  void actions_before_implicit_timestep();

  /// Save the boundary data to file
  void save_boundaries_to_file(ofstream& output_stream, string filename);

  /// Save the solution data to file
  void save_solution_to_file(ofstream& output_stream,
                             string filename,
                             unsigned n_points);

  /// Save the distance data to file
  void save_distance_to_file(ofstream& output_stream, string filename);

  /// Pointer to the "bulk" mesh
  SimpleRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;
  Mesh* ODE_mesh_pt;

  Data* Flux_data_pt;
  Data* Frame_speed_data_pt;
  Data* Frame_travel_data_pt;
};

template<class ELEMENT>
HeleShawChannelProblem<ELEMENT>::HeleShawChannelProblem()
{
  cout << "Problem constructor" << endl;

  this->add_time_stepper_pt(new BDF<2>);

  this->generate_data();

  this->generate_mesh();

  this->assign_mesh();

  this->pin_dirichlet_boundaries();

  this->setup_elements();

  // Setup equation numbering scheme
  cout << "Assign equation numbers." << endl;
  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;
}

/// Set the initial conditions including previous history
template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::set_initial_condition(
  const double& t_step, DocInfo& doc_info)
{
  cout << "Set initial conditions" << endl;

  cout << "Assign initial values" << endl;
  this->assign_initial_values_impulsive(t_step);

  cout << "Set distance history" << endl;
  unsigned element_index = 0;
  unsigned data_index = 0;
  unsigned value_index = 0;

  unsigned n_tstorage = this->ODE_mesh_pt->element_pt(element_index)
                          ->internal_data_pt(data_index)
                          ->time_stepper_pt()
                          ->ntstorage();
  cout << "Storage: " << n_tstorage << endl;
  for (unsigned time_level = 0; time_level < n_tstorage; time_level++)
  {
    cout << time_level << endl;
    double t = this->time_pt()->time(time_level);

    cout << t << endl;
    Vector<double> u0 =
      dynamic_cast<ODEElement*>(this->ODE_mesh_pt->element_pt(element_index))
        ->exact_solution(t);

    cout << u0[0] << endl;
    this->ODE_mesh_pt->element_pt(element_index)
      ->internal_data_pt(data_index)
      ->set_value(time_level, value_index, u0[0]);
  }

  // Solve for initial condition
  this->ODE_mesh_pt->element_pt(element_index)
    ->internal_data_pt(data_index)
    ->pin(value_index);

  this->newton_solve();

  this->ODE_mesh_pt->element_pt(element_index)
    ->internal_data_pt(data_index)
    ->unpin(value_index);

  cout << "Doc solution" << endl;
  this->doc_solution(doc_info);
}

/// Iterate forward in time in steps of t_step until t_final
template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::iterate_timestepper(const double& t_step,
                                                          const double& t_final,
                                                          DocInfo& doc_info)
{
  double t_duration = t_final - this->time_stepper_pt()->time();
  unsigned n_timestep = ceil(t_duration / t_step);
  for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++)
  {
    this->unsteady_newton_solve(t_step);
    this->doc_solution(doc_info);
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  cout << "Doc solution" << endl;

  string data_directory = doc_info.directory();
  ofstream output_stream;
  string filename = "";

  filename = data_directory + "/boundaries.dat";
  this->save_boundaries_to_file(output_stream, filename);

  unsigned n_points = 5;
  filename = data_directory + "/soln" + to_string(doc_info.number()++) + ".dat";
  this->save_solution_to_file(output_stream, filename, n_points);

  filename = data_directory + "/distance.dat";
  this->save_distance_to_file(output_stream, filename);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_data()
{
  double initial_flux = 1.0;
  this->Flux_data_pt = new Data(1);
  this->Flux_data_pt->set_value(0, initial_flux);
  this->Flux_data_pt->pin(0);
  problem_parameter::global_flux_pt = this->Flux_data_pt->value_pt(0);

  double initial_frame_speed = 0.2;
  this->Frame_speed_data_pt = new Data(1);
  this->Frame_speed_data_pt->set_value(0, initial_frame_speed);
  this->Frame_speed_data_pt->pin(0);
  problem_parameter::global_frame_speed_pt =
    this->Frame_speed_data_pt->value_pt(0);

  double initial_frame_travel = 0.0;
  this->Frame_travel_data_pt = new Data(1);
  this->Frame_travel_data_pt->set_value(0, initial_frame_travel);
  this->Frame_travel_data_pt->unpin(0);
  problem_parameter::global_frame_travel_pt =
    this->Frame_travel_data_pt->value_pt(0);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_mesh()
{
  cout << "Generate mesh" << endl;
  unsigned n_x = 20;
  unsigned n_y = 20;
  double l_x = 2.0;
  double l_y = 1.0;

  this->Bulk_mesh_pt = new SimpleRectangularQuadMesh<ELEMENT>(
    n_x, n_y, l_x, l_y, this->time_stepper_pt());

  this->Surface_mesh_pt = new Mesh;

  cout << "Create flux elements" << endl;
  const unsigned flux_boundary = 3;
  this->create_flux_elements(flux_boundary);

  this->ODE_mesh_pt = new Mesh;
  this->ODE_mesh_pt->add_element_pt(
    new ODEElement(this->time_stepper_pt(), new problem_parameter::ODEFunctor));
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::create_flux_elements(
  const unsigned& boundary)
{
  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(boundary);
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(boundary, n));

    int face_index = this->Bulk_mesh_pt->face_index_at_boundary(boundary, n);

    HeleShawFluxElement<ELEMENT>* flux_element_pt =
      new HeleShawFluxElement<ELEMENT>(bulk_element_pt, face_index);

    this->Surface_mesh_pt->add_element_pt(flux_element_pt);
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::assign_mesh()
{
  cout << "Assign mesh" << endl;
  this->add_sub_mesh(this->Bulk_mesh_pt);
  this->add_sub_mesh(this->Surface_mesh_pt);
  this->add_sub_mesh(this->ODE_mesh_pt);

  this->build_global_mesh();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::pin_dirichlet_boundaries()
{
  cout << "Pin Dirichlet boundaries" << endl;
  unsigned n_boundary = this->Bulk_mesh_pt->nboundary();
  bool pin_boundary[n_boundary] = {false};
  pin_boundary[1] = true;
  for (unsigned b = 0; b < n_boundary; b++)
  {
    if (pin_boundary[b])
    {
      cout << "Pinning boundary: " << b << endl;
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_node; n++)
      {
        this->Bulk_mesh_pt->boundary_node_pt(b, n)->pin(0);
      }
    }
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::setup_elements()
{
  cout << "Setup elements" << endl;

  // Find number of elements in mesh
  unsigned n_element = this->Bulk_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->Bulk_mesh_pt->element_pt(i));

    el_pt->upper_wall_fct_pt() = problem_parameter::upper_wall_fct;
  }

  // Find number of elements in mesh
  n_element = this->Surface_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElement<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElement<ELEMENT>*>(
        this->Surface_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = &problem_parameter::get_neumann_bc;
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::set_boundary_conditions()
{
  unsigned n_boundary = this->Bulk_mesh_pt->nboundary();
  bool set_boundary[n_boundary] = {false};
  set_boundary[1] = true;
  for (unsigned b = 0; b < n_boundary; b++)
  {
    if (set_boundary[b])
    {
      cout << "Setting boundary: " << b << endl;
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_node; n++)
      {
        Node* node_pt = this->Bulk_mesh_pt->boundary_node_pt(b, n);
        double value = 0.0;
        Vector<double> x(2);
        x[0] = node_pt->x(0);
        x[1] = node_pt->x(1);
        problem_parameter::get_dirichlet_bc(x, value);
        node_pt->set_value(0, value);
      }
    }
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::actions_before_newton_solve()
{
  this->set_boundary_conditions();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::actions_before_implicit_timestep()
{
  *problem_parameter::global_flux_pt = this->time_pt()->time();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_boundaries_to_file(
  ofstream& output_stream, string filename)
{
  output_stream.open(filename.c_str());
  this->Bulk_mesh_pt->output_boundaries(output_stream);
  output_stream.close();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_solution_to_file(
  ofstream& output_stream, string filename, unsigned n_points)
{
  output_stream.open(filename);
  this->Bulk_mesh_pt->output(output_stream, n_points);
  output_stream.close();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_distance_to_file(
  ofstream& output_stream, string filename)
{
  /// Create and open the file stream
  output_stream.open(filename.c_str(), ofstream::out | ofstream::app);

  /// Write the current time
  unsigned element_index = 0;
  unsigned data_index = 0;
  double t = this->time_pt()->time();
  output_stream << "t: " << t << ", ";

  /// Write the current value
  unsigned value_index = 0;
  double x = this->ODE_mesh_pt->element_pt(element_index)
               ->internal_data_pt(data_index)
               ->value(value_index);
  output_stream << "x: " << x << endl;
  output_stream.close();
}

int main(int argc, char* argv[])
{
  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Label for output
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("data");

  HeleShawChannelProblem<QHeleShawElement<3>> problem;

  cout << "\n\n\nProblem self-test ";
  if (problem.self_test() == 0)
  {
    cout << "passed: Problem can be solved." << std::endl;
  }
  else
  {
    throw OomphLibError(
      "Self test failed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  const double t_step = 0.1;
  const double t_final = 2;

  /// Solve for initial conditions
  problem.set_initial_condition(t_step, doc_info);

  /// Iterate the timestepper once
  problem.iterate_timestepper(t_step, t_final, doc_info);

  return 0;
}
