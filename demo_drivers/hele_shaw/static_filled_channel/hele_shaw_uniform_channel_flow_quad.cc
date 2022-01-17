#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"

using namespace oomph;
using namespace std;

namespace problem_parameter
{
  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    b = 1;
    dbdt = 0;
  }

  void get_dirichlet_bc(const Vector<double>& x, double& p)
  {
    /// At the outlet we set the pressure to be zero
    p = 0.0;
  }

  void get_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the local flux which is dependent on the
    /// upper wall function, inlet_area and total flux
    double total_flux = 1.0;
    double inlet_area = 1.0;
    double dpdx = total_flux / inlet_area;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b);
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
  void create_mesh();

  void create_bulk_mesh();

  void create_surface_mesh();

  /// Create flux elements
  void create_flux_elements(const unsigned& boundary);

  /// Pin dirichlet outlet boundary
  void pin_dirichlet_boundaries();

  /// Upcast elements and finalise setup
  void setup_elements();

  /// Set boundary condition values
  void set_boundary_conditions();

  /// Update the problem specs before solve (empty)
  void actions_before_newton_solve();

  /// Update the problem specs before solve (empty)
  void actions_after_newton_solve() {}

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

  this->create_mesh();

  this->setup_elements();

  this->pin_dirichlet_boundaries();

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
  this->save_solution_to_file(
    output_stream, data_directory + "/soln.dat", n_points);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::create_mesh()
{
  this->create_bulk_mesh();

  this->create_surface_mesh();

  this->add_sub_mesh(this->Bulk_mesh_pt);
  this->add_sub_mesh(this->Surface_mesh_pt);

  this->build_global_mesh();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::create_bulk_mesh()
{
  unsigned n_x = 4;
  unsigned n_y = 4;
  double l_x = 2.0;
  double l_y = 1.0;

  this->Bulk_mesh_pt =
    new SimpleRectangularQuadMesh<ELEMENT>(n_x, n_y, l_x, l_y);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::create_surface_mesh()
{
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

    // Set the flux function pointer
    el_pt->flux_fct_pt() = problem_parameter::get_flux_bc;
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
  cout << "Hele-Shaw channel flow" << endl;

  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Label for output
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("RESLT");

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

  problem.newton_solve();

  problem.doc_solution(doc_info);

  cout << "End of Hele-Shaw channel flow" << endl;

  return 0;
}
