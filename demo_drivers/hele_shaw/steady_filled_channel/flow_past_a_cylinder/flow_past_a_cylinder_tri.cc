#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"

using namespace oomph;
using namespace std;

namespace problem_parameter
{
  double* inlet_integral_pt = 0;
  double* outlet_integral_pt = 0;

  double xcenter = 1;
  double ycenter = 0.5;
  double Major_Radius = 0.2;
  double Minor_Radius = 0.2;
  unsigned circpts = 16;

  /// This is non-dimensionalised to 1
  const double total_flux = 1.0;

  void upper_wall_fct(const Vector<double>& x, double& b, double& dbdt)
  {
    // double tape_height = 0.1;
    // double tape_width = 0.4;
    // double tape_sharpness = 40;
    // double tape_centre_y = 0.5;

    // double y = x[1];

    // b = 1 - tape_height * 0.5 *
    //          (tanh(tape_sharpness * (y - tape_centre_y + 0.5 * tape_width)) -
    //           tanh(tape_sharpness * (y - tape_centre_y - 0.5 * tape_width)));

    double height = 0.0;
    double rms_width = 0.1;
    double centre_x = 0.2;
    double centre_y = 0.5;

    // Transform y such that the domain is between 0 and 1 rather than -1 and 1
    double local_x = x[0];
    double local_y = x[1];
    b = 1.0 - height * std::exp(-(local_x - centre_x) * (local_x - centre_x) /
                                  (2 * rms_width * rms_width) -
                                (local_y - centre_y) * (local_y - centre_y) /
                                  (2 * rms_width * rms_width));
    dbdt = 0.0;
  }

  void get_dirichlet_bc(const Vector<double>& x, double& p)
  {
    /// At the outlet we set the pressure to be zero
    p = 0.0;
  }

  void get_inlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_integral and total flux

    double dpdx = total_flux / *inlet_integral_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b) / 12.0;
  }

  void get_outlet_flux_bc(const Vector<double>& x, double& flux)
  {
    /// At the inlet we set the pressure gradient which is dependent on the
    /// upper wall function, inlet_integral and total flux

    double dpdx = total_flux / *outlet_integral_pt;

    double b;
    double dbdt;
    upper_wall_fct(x, b, dbdt);

    flux = dpdx * (b * b * b) / 12.0;
  }

  enum
  {
    LOWER_BOUNDARY,
    OUTLET_BOUNDARY,
    UPPER_BOUNDARY,
    INLET_BOUNDARY,
    FIRST_BUBBLE_BOUNDARY,
    SECOND_BUBBLE_BOUNDARY
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

  /// Doc the solution
  void doc_solution(DocInfo& doc_info);

private:
  /// Generate mesh
  void generate_mesh();

  /// Generate bulk mesh
  void generate_bulk_mesh();

  /// Generate rectangular outer mesh boundaries
  void generate_rect_boundary();

  /// Generate ellipse inner mesh boundary
  void generate_inner_boundary();

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

  /// Pointer to the outer boundary polyline
  TriangleMeshPolygon* Rect_boundary_polyline_pt;

  /// Pointer to the inner boundary polyline
  TriangleMeshPolygon* Cylinder_boundary_polyline_pt;

  /// Pointer to the "bulk" mesh
  // TriangleMesh<ELEMENT>* Bulk_mesh_pt;
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
  this->generate_rect_boundary();

  this->generate_inner_boundary();

  // Convert to "closed curve" objects
  TriangleMeshClosedCurve* rect_closed_curve_pt =
    this->Rect_boundary_polyline_pt;

  Vector<TriangleMeshClosedCurve*> cylinder_closed_curve_pt(1);
  cylinder_closed_curve_pt[0] = this->Cylinder_boundary_polyline_pt;

  // Generate mesh parameters for external mesh generator "Triangle"
  TriangleMeshParameters triangle_mesh_parameters(rect_closed_curve_pt);
  triangle_mesh_parameters.internal_closed_curve_pt() = cylinder_closed_curve_pt;

  double maximum_default_element_area = 1e-2;
  triangle_mesh_parameters.element_area() = maximum_default_element_area;

  // Call external mesh generator
  this->Bulk_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_rect_boundary()
{
  double l_x = 2.0;
  double l_y = 1.0;

  cout << "Generate boundary" << endl;
  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

  cout << "Create 2 vertex vector" << endl;
  Vector<Vector<double>> vertex_coord(2);
  for (unsigned i = 0; i < 2; i++)
  {
    vertex_coord[i].resize(2);
  }

  /// First vertex
  vertex_coord[0][0] = 0;
  vertex_coord[0][1] = 0;
  /// Second vertex
  vertex_coord[1][0] = l_x;
  vertex_coord[1][1] = 0;

  cout << "Lower boundary, ";
  boundary_polyline_pt[0] =
    new TriangleMeshPolyLine(vertex_coord, problem_parameter::LOWER_BOUNDARY);

  /// First vertex
  vertex_coord[0][0] = l_x;
  vertex_coord[0][1] = 0;
  /// Second vertex
  vertex_coord[1][0] = l_x;
  vertex_coord[1][1] = l_y;

  cout << "outer boundary, ";
  boundary_polyline_pt[1] =
    new TriangleMeshPolyLine(vertex_coord, problem_parameter::OUTLET_BOUNDARY);

  /// First vertex
  vertex_coord[0][0] = l_x;
  vertex_coord[0][1] = l_y;
  /// Second vertex
  vertex_coord[1][0] = 0;
  vertex_coord[1][1] = l_y;

  cout << "upper boundary, ";
  boundary_polyline_pt[2] =
    new TriangleMeshPolyLine(vertex_coord, problem_parameter::UPPER_BOUNDARY);

  /// First vertex
  vertex_coord[0][0] = 0;
  vertex_coord[0][1] = l_y;
  /// Second vertex
  vertex_coord[1][0] = 0;
  vertex_coord[1][1] = 0;

  cout << "inlet boundary, " << endl;
  boundary_polyline_pt[3] =
    new TriangleMeshPolyLine(vertex_coord, problem_parameter::INLET_BOUNDARY);

  cout << "Create triangle mesh polygon" << endl;
  // Create the triangle mesh polygon for rectangle boundary
  Rect_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_inner_boundary()
{
  // Place it smack in the middle of the channel
  double x_center = problem_parameter::xcenter;
  double y_center = problem_parameter::ycenter;
  double major_radius = problem_parameter::Major_Radius;
  double minor_radius = problem_parameter::Minor_Radius;
  Ellipse* bubble_pt = new Ellipse(major_radius, minor_radius);

  // Intrinsic coordinate along GeomObject defining the bubble
  Vector<double> zeta(1);

  // Position vector to GeomObject defining the bubble
  Vector<double> coord(2);

  // Number of points defining bubble
  unsigned npoints = problem_parameter::circpts; // 16;
  double unit_zeta = MathematicalConstants::Pi / double(npoints - 1);

  // This bubble is bounded by two distinct boundaries, each
  // represented by its own polyline
  Vector<TriangleMeshCurveSection*> bubble_polyline_pt(2);

  // Vertex coordinates
  Vector<Vector<double>> bubble_vertex(npoints);

  // Create points on boundary
  for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
  {
    // Resize the vector
    bubble_vertex[ipoint].resize(2);

    // Get the coordinates
    zeta[0] = unit_zeta * double(ipoint);
    bubble_pt->position(zeta, coord);

    // Shift
    bubble_vertex[ipoint][0] = coord[0] + x_center;
    bubble_vertex[ipoint][1] = coord[1] + y_center;
  }

  // Build the 1st bubble polyline
  bubble_polyline_pt[0] = new TriangleMeshPolyLine(
    bubble_vertex, problem_parameter::FIRST_BUBBLE_BOUNDARY);

  // Second boundary of bubble
  for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
  {
    // Resize the vector
    bubble_vertex[ipoint].resize(2);

    // Get the coordinates
    zeta[0] = (unit_zeta * double(ipoint)) + MathematicalConstants::Pi;
    bubble_pt->position(zeta, coord);

    // Shift
    bubble_vertex[ipoint][0] = coord[0] + x_center;
    bubble_vertex[ipoint][1] = coord[1] + y_center;
  }

  // Build the 2nd bubble polyline
  bubble_polyline_pt[1] = new TriangleMeshPolyLine(
    bubble_vertex, problem_parameter::SECOND_BUBBLE_BOUNDARY);


  // Define coordinates of a point inside the bubble
  Vector<double> bubble_center(2);
  bubble_center[0] = x_center;
  bubble_center[1] = y_center;

  bubble_polyline_pt[0]->set_maximum_length(0.02);
  bubble_polyline_pt[1]->set_maximum_length(0.02);

  // Create closed polygon from two polylines
  Cylinder_boundary_polyline_pt =
    new TriangleMeshPolygon(bubble_polyline_pt, bubble_center);
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

    HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
      new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
        bulk_element_pt,
        face_index,
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

    HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
      new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
        bulk_element_pt,
        face_index,
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
    HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
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
    HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
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
  HeleShawChannelProblem<THeleShawElement<3>> problem;

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
