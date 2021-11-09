#include <cmath>
#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"

using namespace oomph;
using namespace std;

namespace problem_parameter
{
  double* global_time_pt = 0;
  const double speed = 1.0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    double tape_height = 0.1;
    double tape_width = 0.4;
    double tape_sharpness = 40;
    double tape_centre_y = 0.5;
    double tape_start_x = 2.5 - (*global_time_pt) * speed;

    double y = x[1];

    b =
      1.0 - tape_height * 0.5 *
              (tanh(tape_sharpness * (y - tape_centre_y + 0.5 * tape_width)) -
               tanh(tape_sharpness * (y - tape_centre_y - 0.5 * tape_width))) *
              0.5 * (tanh(tape_sharpness * (x[0] - tape_start_x)) + 1);

    double sech2 = 1 - tanh(tape_sharpness * (x[0] - tape_start_x)) *
                         tanh(tape_sharpness * (x[0] - tape_start_x));

    dbdt = tape_height * 0.5 *
           (tanh(tape_sharpness * (y - tape_centre_y + 0.5 * tape_width)) -
            tanh(tape_sharpness * (y - tape_centre_y - 0.5 * tape_width))) *
           0.5 * speed * tape_sharpness * sech2;
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
    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);
    double total_flux = 1.0;
    /// This needs to be calculated
    double inlet_area = 1.0;
    dpdx = total_flux / inlet_area * (b * b * b);
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

  /// Pointer to the "bulk" mesh
  SimpleRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;
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

  this->Surface_mesh_pt = new Mesh;

  cout << "Create flux elements" << endl;
  const unsigned flux_boundary = 3;
  this->create_flux_elements(flux_boundary);
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
  double tF = 1;
  /// Iterate the timestepper using the fixed time step until the final time
  problem.iterate_timestepper(dt, tF, doc_info);

  return 0;
}
