#include <cmath>
#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"

using namespace oomph;
using namespace std;

namespace problem_parameter
{
  double* global_time_pt = 0;
  double* inlet_b3_pt = 0;
  double* outlet_b3_pt = 0;
  const double total_flux = 1.0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    double perturbation_height = 0.5;
    double rms_width = 0.1;
    double centre_x =
      0.5 - 3.0 * (pow(*global_time_pt, 2.0) / 2.0 - pow(*global_time_pt, 3.0) / 3.0);
    double dcentre_xdt = -3.0 * *global_time_pt * (1.0 - *global_time_pt);
    double centre_y = 0.5;

    // Transform y such that the domain is between 0 and 1 rather than -1 and 1
    double local_x = x[0];
    double local_y = x[1];

    double f = 0.0;
    f = -(local_x - centre_x) * (local_x - centre_x) /
          (2.0 * rms_width * rms_width) -
        (local_y - centre_y) * (local_y - centre_y) /
          (2.0 * rms_width * rms_width);
    double dfdt = 0.0;
    dfdt =
      -2.0 * dcentre_xdt * (local_x - centre_x) / (2.0 * rms_width * rms_width);

    b = 1.0 - perturbation_height * exp(f);

    dbdt = -dfdt * perturbation_height * exp(f);
  }

  void get_dirichlet_bc(const Vector<double>& x, double& p)
  {
    p = 0.0;
  }

  void get_inlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_b3 and total flux

    double G = total_flux / *inlet_b3_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = -G * pow(b, 3.0) / 12;
  }

  void get_outlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_b3 and total flux

    double G = total_flux / *outlet_b3_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = G * pow(b, 3.0) / 12;
  }


} // namespace problem_parameter

template<class ELEMENT>
class HeleShawChannelProblem : public Problem
{
public:
  /// Constructor
  HeleShawChannelProblem();

  /// Destructor (empty)
  ~HeleShawChannelProblem() {}

  /// Use a newton solve to set the initial conditions
  void solve_for_initial_conditions(DocInfo& doc_info);

  /// Iterate forward in time
  void iterate_timestepper(const double& t_step,
                           const double& t_final,
                           DocInfo& doc_info);

  /// Doc the solution
  void doc_solution(DocInfo& doc_info);

private:
  /// Generate mesh
  void generate_mesh();

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

  /// Pointer to the bulk mesh
  SimpleRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the inlet meshes
  Mesh* Inlet_surface_mesh_pt;
  Mesh* Outlet_surface_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;

  /// Pointers to the inlet and outlet integral data
  Data* Inlet_integral_data_pt;
  Data* Outlet_integral_data_pt;

  enum
  {
    Lower_boundary,
    Outlet_boundary,
    Upper_boundary,
    Inlet_boundary
  };
};

template<class ELEMENT>
HeleShawChannelProblem<ELEMENT>::HeleShawChannelProblem()
{
  cout << "Problem constructor" << endl;

  this->add_time_stepper_pt(new BDF<1>);
  problem_parameter::global_time_pt = &this->time_pt()->time();

  this->generate_mesh();

  this->assign_mesh();

  this->pin_dirichlet_boundaries();

  this->setup_elements();

  // Setup equation numbering scheme
  cout << "Assign equation numbers." << endl;
  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  cout << "Doc solution" << endl;

  string data_directory = doc_info.directory();
  ofstream output_stream;

  string filename = data_directory + "/boundaries.dat";
  this->save_boundaries_to_file(output_stream, filename);

  unsigned n_points = 5;
  this->save_solution_to_file(output_stream,
                              data_directory + "/soln" +
                                to_string(doc_info.number()) + ".dat",
                              n_points);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_mesh()
{
  cout << "Generate mesh" << endl;
  unsigned n_x = 32;
  unsigned n_y = 32;
  double l_x = 2.0;
  double l_y = 1.0;

  this->Bulk_mesh_pt = new SimpleRectangularQuadMesh<ELEMENT>(
    n_x, n_y, l_x, l_y, this->time_stepper_pt());

  cout << "create inlet data" << endl;
  /// Pointer to inlet integral data
  unsigned number_of_values = 1;
  unsigned index = 0;
  Inlet_integral_data_pt = new Data(number_of_values);
  Outlet_integral_data_pt = new Data(number_of_values);

  Inlet_integral_data_pt->set_value(index, 1.0);
  Inlet_integral_data_pt->unpin(index);
  problem_parameter::inlet_b3_pt = Inlet_integral_data_pt->value_pt(index);

  Outlet_integral_data_pt->set_value(index, 1.0);
  Outlet_integral_data_pt->unpin(index);
  problem_parameter::outlet_b3_pt = Outlet_integral_data_pt->value_pt(index);

  cout << "Add to info mesh" << endl;
  this->Info_mesh_pt = new Mesh;
  this->Info_mesh_pt->add_element_pt(new InfoElement(Inlet_integral_data_pt));
  this->Info_mesh_pt->add_element_pt(new InfoElement(Outlet_integral_data_pt));

  cout << "Create flux elements" << endl;
  this->Inlet_surface_mesh_pt = new Mesh;
  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(Inlet_boundary);
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(Inlet_boundary, n));

    int face_index =
      this->Bulk_mesh_pt->face_index_at_boundary(Inlet_boundary, n);

    HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
      new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
        bulk_element_pt, face_index, Inlet_integral_data_pt);

    this->Inlet_surface_mesh_pt->add_element_pt(flux_element_pt);
  }

  this->Outlet_surface_mesh_pt = new Mesh;
  n_element = this->Bulk_mesh_pt->nboundary_element(Outlet_boundary);
  cout << n_element << endl;
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(Outlet_boundary, n));

    int face_index =
      this->Bulk_mesh_pt->face_index_at_boundary(Outlet_boundary, n);

    HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
      new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
        bulk_element_pt, face_index, Outlet_integral_data_pt);

    this->Outlet_surface_mesh_pt->add_element_pt(flux_element_pt);
  }
  cout << "Done" << endl;
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::assign_mesh()
{
  cout << "Assign mesh" << endl;
  this->add_sub_mesh(this->Bulk_mesh_pt);
  this->add_sub_mesh(this->Inlet_surface_mesh_pt);
  this->add_sub_mesh(this->Outlet_surface_mesh_pt);
  this->add_sub_mesh(this->Info_mesh_pt);

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
  n_element = this->Inlet_surface_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
        this->Inlet_surface_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = problem_parameter::get_inlet_flux_bc;
  }

  // Find number of elements in mesh
  n_element = this->Outlet_surface_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
        this->Outlet_surface_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = problem_parameter::get_inlet_flux_bc;
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
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::solve_for_initial_conditions(
  DocInfo& doc_info)
{
  this->newton_solve();
  this->doc_solution(doc_info);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::iterate_timestepper(const double& t_step,
                                                          const double& t_final,
                                                          DocInfo& doc_info)
{
  unsigned n_timestep = ceil(t_final / t_step);

  for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++)
  {
    cout << "t: " << this->time_pt()->time() << endl;

    this->unsteady_newton_solve(t_step);

    doc_info.number()++;
    this->doc_solution(doc_info);
  }
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
  unsigned index = 0;
  double my_integral =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value(index);
  output_stream << "Inflow Integral = " << my_integral << endl;
  output_stream.close();
}

int main(int argc, char* argv[])
{
  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Label for output
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT/");

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

  /// Solve for initial conditions
  problem.solve_for_initial_conditions(doc_info);

  double dt = 0.1;
  double tF = 1.0;

  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  return 0;
}
