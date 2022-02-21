#ifndef OOMPH_FILLED_CHANNEL_FLOW_PROBLEM_HEADER
#define OOMPH_FILLED_CHANNEL_FLOW_PROBLEM_HEADER

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"
#include "local_pert_parameter.h"
//#include "tape_parameter.h"
//#include "problem_parameter.h"

template<class ELEMENT>
class HeleShawChannelProblem : public Problem
{
public:
  /// Constructor
  HeleShawChannelProblem() {}

  /// Destructor (empty)
  ~HeleShawChannelProblem() {}


  void setup();

  /// Doc the solution
  void doc_solution(DocInfo& doc_info);

protected:
  /// Pointer to the "bulk" mesh
  Mesh* Bulk_mesh_pt;

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

  /// Pointer to the "surface" mesh
  Mesh* Inlet_surface_mesh_pt;
  Mesh* Outlet_surface_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;
};

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::setup()
{
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
  generate_bulk_mesh();
  cout << "finished generate_bulk_mesh" << endl;
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
  double domain_length = problem_parameter::domain_length;
  double domain_width = problem_parameter::domain_width;

  unsigned id = 0;
  /// Create a mesh curve section for each boundary
  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

  // Each polyline only has two vertices -- provide storage for their
  // coordinates
  Vector<Vector<double>> vertex_coord(2);
  for (unsigned i = 0; i < 2; i++)
  {
    vertex_coord[i].resize(2);
  }

  // First polyline: bottom
  vertex_coord[0][0] = 0.0;
  vertex_coord[0][1] = 0.0;
  vertex_coord[1][0] = domain_length;
  vertex_coord[1][1] = 0.0;

  // Build the 1st boundary polyline
  boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord, id);
  id++;

  // Second boundary polyline: right
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = domain_length;
  vertex_coord[1][1] = domain_width;

  // Build the 2nd boundary polyline
  boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord, id);
  id++;

  // Third boundary polyline: top
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = 0.0;
  vertex_coord[1][1] = domain_width;

  // Build the 3rd boundary polyline
  boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord, id);
  id++;

  // Fourth boundary polyline: left
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = 0.0;
  vertex_coord[1][1] = 0.0;

  // Build the 4th boundary polyline
  boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord, id);
  id++;

  // Create the triangle mesh polygon for outer boundary
  TriangleMeshPolygon* Outer_boundary_polygon_pt =
    new TriangleMeshPolygon(boundary_polyline_pt);

  // Convert to "closed curve" objects
  TriangleMeshClosedCurve* rect_closed_curve_pt = Outer_boundary_polygon_pt;

  // Generate mesh parameters for external mesh generator "Triangle"
  TriangleMeshParameters triangle_mesh_parameters(rect_closed_curve_pt);

  double maximum_default_element_area = 1e-2;
  triangle_mesh_parameters.element_area() = maximum_default_element_area;

  // Call external mesh generator
  this->Bulk_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);
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
  problem_parameter::inlet_integral_pt =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value_pt(index);
  problem_parameter::outlet_integral_pt =
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
  for (unsigned index = 0; index < 1; index++)
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

#endif
