#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"

using namespace oomph;
using namespace std;

namespace problem_parameter
{
  double* inlet_area_pt = 0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    double tape_height = 0.1;
    double tape_width = 0.4;
    double tape_sharpness = 40;
    double tape_centre_y = 0.5;

    double y = x[1];

    b = 1 - tape_height * 0.5 *
              (tanh(tape_sharpness * (y - tape_centre_y + 0.5 * tape_width)) -
               tanh(tape_sharpness * (y - tape_centre_y - 0.5 * tape_width)));
    dbdt = 0.0;
  }

  void get_dirichlet_bc(const Vector<double>& x, double& p)
  {
    /// At the outlet we set the pressure to be zero
    p = 0.0;
  }

  void get_neumann_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_area and total flux

    /// This is non-dimensionalised to 1
    double total_flux = 1.0;
    double dpdx = total_flux / *inlet_area_pt;

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
  void generate_mesh();

  /// Generate bulk mesh
  void generate_bulk_mesh();

  /// Generate info mesh
  void generate_info_mesh();

  /// Generate surface elements
  void generate_surface_mesh();

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
  SimpleRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;
};

template<class ELEMENT>
HeleShawChannelProblem<ELEMENT>::HeleShawChannelProblem()
{
  this->generate_mesh();

  this->upcast_and_finalise_elements();

  this->pin_data();

  this->set_boundary_conditions();

  // Setup equation numbering scheme
  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  string data_directory = doc_info.directory();
  ofstream output_stream;

  string filename = data_directory + "/boundaries.dat";
  this->save_boundaries_to_file(output_stream, filename);

  unsigned n_points = 5;
  this->save_solution_to_file(
    output_stream, data_directory + "/soln.dat", n_points);

  filename = data_directory + "/integral.dat";
  this->save_integral_to_file(output_stream, filename);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_mesh()
{
  this->generate_bulk_mesh();
  this->generate_info_mesh();
  this->generate_surface_mesh();

  this->add_sub_mesh(this->Bulk_mesh_pt);
  this->add_sub_mesh(this->Surface_mesh_pt);
  this->add_sub_mesh(this->Info_mesh_pt);

  this->build_global_mesh();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_bulk_mesh()
{
  unsigned n_x = 4;
  unsigned n_y = 4;
  double l_x = 2.0;
  double l_y = 1.0;

  this->Bulk_mesh_pt =
    new SimpleRectangularQuadMesh<ELEMENT>(n_x, n_y, l_x, l_y);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_info_mesh()
{
  unsigned number_of_values = 1;

  this->Info_mesh_pt = new Mesh;

  Data* Inlet_integral_data_pt = new Data(number_of_values);
  this->Info_mesh_pt->add_element_pt(new InfoElement(Inlet_integral_data_pt));
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_surface_mesh()
{
  const unsigned boundary = 3;

  this->Surface_mesh_pt = new Mesh;

  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(boundary);
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(boundary, n));

    int face_index = this->Bulk_mesh_pt->face_index_at_boundary(boundary, n);

    HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
      new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
        bulk_element_pt,
        face_index,
        this->Info_mesh_pt->element_pt(0)->internal_data_pt(0));

    this->Surface_mesh_pt->add_element_pt(flux_element_pt);
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
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::upcast_and_finalise_elements()
{
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

  unsigned index = 0;
  problem_parameter::inlet_area_pt =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value_pt(index);

  // Find number of elements in mesh
  n_element = this->Surface_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
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

  unsigned index = 0;
  /// Integral must be initialised to something non-zero
  this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->set_value(index, 1.0);
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
  unsigned index = 0;
  double my_integral =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value(index);
  output_stream << "Integral = " << my_integral << endl;
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

  /// Call problem solve
  problem.newton_solve();

  /// Document solution
  problem.doc_solution(doc_info);

  return 0;
}
