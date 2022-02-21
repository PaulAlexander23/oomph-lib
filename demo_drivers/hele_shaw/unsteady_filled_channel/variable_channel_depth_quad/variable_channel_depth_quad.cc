#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"

using namespace oomph;
using namespace std;

namespace problem_parameter
{
  /// Global data
  double* global_time_pt = 0;
  double* inlet_b3_pt = 0;
  double* outlet_b3_pt = 0;

  /// This is non-dimensionalised to 1
  const double total_flux = 1.0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    double perturbation_height = 0.5;
    double rms_width = 0.1;

    double period = 1.1;
    double centre_x = 1 + 0.5 * sin(8 * atan(1) * (*global_time_pt) / period);
    double dcentre_xdt = 0.5 * 8 * atan(1) / period *
                         cos(8 * atan(1) * (*global_time_pt) / period);
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
    /// At the outlet we set the pressure to be zero
    p = 0.0;
  }

  void get_inlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_area and total flux
    double dpdx = total_flux / *inlet_b3_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b) / 12.0;
  }

  void get_outlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_area and total flux
    double dpdx = total_flux / *outlet_b3_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b) / 12.0;
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

  /// Doc the solution
  void doc_solution(DocInfo& doc_info);

private:
  /// Generate mesh
  void generate_mesh();

  /// Generate bulk mesh
  void generate_bulk_mesh();

  /// Generate info mesh
  void generate_info_mesh();

  /// Generate surface elements
  void generate_inlet_surface_mesh();
  void generate_outlet_surface_mesh();

  /// Pin dirichlet outlet boundary
  void pin_data();

  /// Upcast elements and finalise setup
  void upcast_and_finalise_elements();

  /// Set boundary condition values
  void set_boundary_conditions();

  /// Save the boundary data to file
  void save_boundaries_to_file(ofstream& output_stream, string filename);

  /// Save the solution data to file
  void save_solution_to_file(ofstream& output_stream,
                             string filename,
                             unsigned n_points);

  /// Save the integral to file
  void save_integral_to_file(ofstream& output_stream, string filename);

  /// Pointer to the "bulk" mesh
  // SimpleRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;
  Mesh* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Inlet_surface_mesh_pt;
  Mesh* Outlet_surface_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;
};

template<class ELEMENT>
HeleShawChannelProblem<ELEMENT>::HeleShawChannelProblem()
{
  this->add_time_stepper_pt(new BDF<1>);
  problem_parameter::global_time_pt = &this->time_pt()->time();

  cout << "generate_mesh" << endl;
  this->generate_mesh();

  cout << "upcast_and_finalise_elements" << endl;
  this->upcast_and_finalise_elements();

  cout << "pin_data" << endl;
  this->pin_data();

  cout << "set_boundary_conditions" << endl;
  this->set_boundary_conditions();

  // Setup equation numbering scheme
  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  string data_directory = doc_info.directory();
  string doc_num = to_string(doc_info.number());
  ofstream output_stream;

  string filename = data_directory + "/boundaries" + doc_num + ".dat";
  this->save_boundaries_to_file(output_stream, filename);

  unsigned n_points = 5;
  this->save_solution_to_file(
    output_stream, data_directory + "/soln" + doc_num + ".dat", n_points);

  filename = data_directory + "/integral" + doc_num + ".dat";
  this->save_integral_to_file(output_stream, filename);

  doc_info.number()++;
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_mesh()
{
  this->generate_bulk_mesh();
  this->generate_info_mesh();
  this->generate_inlet_surface_mesh();
  this->generate_outlet_surface_mesh();

  this->add_sub_mesh(this->Bulk_mesh_pt);
  this->add_sub_mesh(this->Inlet_surface_mesh_pt);
  this->add_sub_mesh(this->Outlet_surface_mesh_pt);
  this->add_sub_mesh(this->Info_mesh_pt);

  this->build_global_mesh();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_bulk_mesh()
{
  unsigned n_x = 20;
  unsigned n_y = 20;
  double l_x = 2.0;
  double l_y = 1.0;

  this->Bulk_mesh_pt =
    new SimpleRectangularQuadMesh<ELEMENT>(n_x, n_y, l_x, l_y);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_info_mesh()
{
  this->Info_mesh_pt = new Mesh;

  unsigned number_of_values = 1;
  Data* Inlet_integral_data_pt = new Data(number_of_values);
  this->Info_mesh_pt->add_element_pt(new InfoElement(Inlet_integral_data_pt));
  Data* Outlet_integral_data_pt = new Data(number_of_values);
  this->Info_mesh_pt->add_element_pt(new InfoElement(Outlet_integral_data_pt));
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_inlet_surface_mesh()
{
  const unsigned boundary = 3;

  this->Inlet_surface_mesh_pt = new Mesh;

  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(boundary);
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(boundary, n));

    int face_index = this->Bulk_mesh_pt->face_index_at_boundary(boundary, n);

    HeleShawFluxElement<ELEMENT>* flux_element_pt =
      new HeleShawFluxElement<ELEMENT>(bulk_element_pt, face_index);

    flux_element_pt->add_b3_data_pt(
      this->Info_mesh_pt->element_pt(0)->internal_data_pt(0));

    this->Inlet_surface_mesh_pt->add_element_pt(flux_element_pt);
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_outlet_surface_mesh()
{
  const unsigned boundary = 1;

  this->Outlet_surface_mesh_pt = new Mesh;

  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(boundary);
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(boundary, n));

    int face_index = this->Bulk_mesh_pt->face_index_at_boundary(boundary, n);

    HeleShawFluxElement<ELEMENT>* flux_element_pt =
      new HeleShawFluxElement<ELEMENT>(bulk_element_pt, face_index);

    flux_element_pt->add_b3_data_pt(
      this->Info_mesh_pt->element_pt(1)->internal_data_pt(0));

    this->Outlet_surface_mesh_pt->add_element_pt(flux_element_pt);
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::pin_data()
{
  unsigned n_boundary = this->Bulk_mesh_pt->nboundary();
  bool pin_boundary[n_boundary] = {false};
  pin_boundary[1] = true;
  for (unsigned b = 0; b < n_boundary; b++)
  {
    if (pin_boundary[b])
    {
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_node; n++)
      {
        this->Bulk_mesh_pt->boundary_node_pt(b, n)->pin(0);
      }
    }
  }

  unsigned index = 0;
  this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->unpin(index);
  this->Info_mesh_pt->element_pt(1)->internal_data_pt(0)->unpin(index);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::upcast_and_finalise_elements()
{
  /// Bulk mesh
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

  /// Info mesh
  unsigned index = 0;
  problem_parameter::inlet_b3_pt =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value_pt(index);
  problem_parameter::outlet_b3_pt =
    this->Info_mesh_pt->element_pt(1)->internal_data_pt(0)->value_pt(index);

  /// Inlet mesh
  // Find number of elements in mesh
  n_element = this->Inlet_surface_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElement<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElement<ELEMENT>*>(
        this->Inlet_surface_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = &problem_parameter::get_inlet_flux_bc;
  }

  /// Outlet mesh
  // Find number of elements in mesh
  n_element = this->Outlet_surface_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElement<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElement<ELEMENT>*>(
        this->Outlet_surface_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = &problem_parameter::get_outlet_flux_bc;
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::set_boundary_conditions()
{
  // unsigned n_boundary = this->Bulk_mesh_pt->nboundary();
  // bool set_boundary[n_boundary] = {false};
  // set_boundary[1] = true;
  // for (unsigned b = 0; b < n_boundary; b++)
  //{
  //  if (set_boundary[b])
  //  {
  //    unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
  //    for (unsigned n = 0; n < n_node; n++)
  //    {
  //      Node* node_pt = this->Bulk_mesh_pt->boundary_node_pt(b, n);
  //      double value = 0.0;
  //      Vector<double> x(2);
  //      x[0] = node_pt->x(0);
  //      x[1] = node_pt->x(1);
  //      problem_parameter::get_dirichlet_bc(x, value);
  //      node_pt->set_value(0, value);
  //    }
  //  }
  //}

  /// Pin a single node in the bulk mesh
  /// Pinning the zeroth node on the zeroth boundary
  const unsigned b = 0;
  const unsigned n = 0;
  Node* node_pt = this->Bulk_mesh_pt->boundary_node_pt(b, n);
  double value = 0.0;
  Vector<double> x(2);
  x[0] = node_pt->x(0);
  x[1] = node_pt->x(1);
  problem_parameter::get_dirichlet_bc(x, value);
  node_pt->set_value(0, value);

  cout << "Integral info mesh" << endl;
  /// Integral info mesh
  const unsigned value_index = 0;
  for (unsigned index = 0; index < 2; index++)
  {
    cout << "index: " << index << endl;
    /// Integral value must be initialised to something non-zero
    this->Info_mesh_pt->element_pt(index)->internal_data_pt(0)->set_value(
      value_index, 1.0);
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

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_integral_to_file(
  ofstream& output_stream, string filename)
{
  output_stream.open(filename.c_str());
  double my_integral = 0.0;
  for (unsigned index = 0; index < 2; index++)
  {
    my_integral =
      this->Info_mesh_pt->element_pt(index)->internal_data_pt(0)->value(0);
    output_stream << "Integral " << index << "= " << my_integral << endl;
  }
  output_stream.close();
}

int main(int argc, char* argv[])
{
  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT");

  /// Create problem object
  HeleShawChannelProblem<QHeleShawElement<3>> problem;

  /// Run problem self test
  if (problem.self_test())
  {
    throw OomphLibError(
      "Self test failed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  // Solve and document the initial conditions
  problem.newton_solve();
  problem.doc_solution(doc_info);

  double time_step = 0.1;
  double time_final = 1.0;

  unsigned n_timestep = ceil(time_final / time_step);
  for (unsigned n = 0; n < n_timestep; n++)
  {
    // Iterate the timestepper once
    problem.unsteady_newton_solve(time_step);
    problem.doc_solution(doc_info);
  }

  return 0;
}
