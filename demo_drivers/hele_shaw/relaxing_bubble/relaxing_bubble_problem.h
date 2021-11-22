#ifndef OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER
#define OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "relaxing_bubble_parameters.h"

template<class ELEMENT>
class RelaxingBubbleProblem : public Problem
{
  /// Variables
public:
  /// Enumeration of mesh boundaries
  enum
  {
    Left_boundary_id = 0,
    Top_boundary_id = 1,
    Right_boundary_id = 2,
    Bottom_boundary_id = 3,
    First_bubble_boundary_id = 4,
    Second_bubble_boundary_id = 5,
  };

private:
  TriangleMeshPolygon* Outer_boundary_polygon_pt;
  TriangleMeshPolygon* Surface_polygon_pt;

  RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;
  Mesh* Surface_mesh_pt;


  /// Functions
public:
  /// Constructor
  RelaxingBubbleProblem();

  /// Destructor
  ~RelaxingBubbleProblem() {}

  /// Use a newton solve to set the initial conditions
  void solve_for_initial_conditions(DocInfo& doc_info) {}

  /// Iterate forward in time
  void iterate_timestepper(const double& t_step,
                           const double& t_final,
                           DocInfo& doc_info)
  {
  }

  /// Doc the solution
  void doc_solution(DocInfo& doc_info) {}

private:
  void generate_mesh();
  void generate_outer_boundary_polygon();
  void generate_surface_polygon();
  void generate_fluid_mesh();
  void generate_surface_mesh();
};

template<class ELEMENT>
RelaxingBubbleProblem<ELEMENT>::RelaxingBubbleProblem()
{
  bool adaptive_timestepping = false;
  add_time_stepper_pt(new BDF<2>(adaptive_timestepping));

  generate_outer_boundary_polygon();
  generate_surface_polygon();

  generate_mesh();
}

template<class ELEMENT>
void RelaxingBubbleProblem<ELEMENT>::generate_outer_boundary_polygon()
{
  cout << "Generate outer boundary polygon" << endl;
  double domain_length = 1.0;
  double domain_width = 1.0;

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
  boundary_polyline_pt[0] =
    new TriangleMeshPolyLine(vertex_coord, Bottom_boundary_id);

  // Second boundary polyline: right
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = domain_length;
  vertex_coord[1][1] = domain_width;

  // Build the 2nd boundary polyline
  boundary_polyline_pt[1] =
    new TriangleMeshPolyLine(vertex_coord, Right_boundary_id);

  // Third boundary polyline: top
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = 0.0;
  vertex_coord[1][1] = domain_width;

  // Build the 3rd boundary polyline
  boundary_polyline_pt[2] =
    new TriangleMeshPolyLine(vertex_coord, Top_boundary_id);

  // Fourth boundary polyline: left
  vertex_coord[0][0] = vertex_coord[1][0];
  vertex_coord[0][1] = vertex_coord[1][1];
  vertex_coord[1][0] = 0.0;
  vertex_coord[1][1] = 0.0;

  // Build the 4th boundary polyline
  boundary_polyline_pt[3] =
    new TriangleMeshPolyLine(vertex_coord, Left_boundary_id);

  // Create the triangle mesh polygon for outer boundary
  Outer_boundary_polygon_pt = new TriangleMeshPolygon(boundary_polyline_pt);
}

template<class ELEMENT>
void RelaxingBubbleProblem<ELEMENT>::generate_surface_polygon()
{
  cout << "Generate surface polygon" << endl;
  double x_center = 0.5;
  double y_center = 0.5;
  double major_radius = 0.4;
  double minor_radius = 0.2;
  unsigned npoints = 16;
  double zeta_step = MathematicalConstants::Pi / double(npoints - 1);

  // Intrinsic coordinate along GeomObject defining the bubble
  Vector<double> zeta(1);

  // Position vector to GeomObject defining the bubble
  Vector<double> coord(2);

  Ellipse* ellipse_pt = new Ellipse(major_radius, minor_radius);

  /// Triangle requires at least two polylines for a single polygon
  Vector<TriangleMeshCurveSection*> surface_polyline_pt(2);

  // Vertex coordinates
  Vector<Vector<double>> bubble_vertex(npoints);

  for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
  {
    bubble_vertex[ipoint].resize(2);

    // Get the coordinates
    zeta[0] = zeta_step * double(ipoint);
    ellipse_pt->position(zeta, coord);

    // Shift
    bubble_vertex[ipoint][0] = coord[0] + x_center;
    bubble_vertex[ipoint][1] = coord[1] + y_center;
  }

  // Build the 1st bubble polyline
  surface_polyline_pt[0] =
    new TriangleMeshPolyLine(bubble_vertex, First_bubble_boundary_id);

  for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
  {
    bubble_vertex[ipoint].resize(2);

    // Get the coordinates
    zeta[0] = zeta_step * double(ipoint) + MathematicalConstants::Pi;
    ellipse_pt->position(zeta, coord);

    // Shift
    bubble_vertex[ipoint][0] = coord[0] + x_center;
    bubble_vertex[ipoint][1] = coord[1] + y_center;
  }

  // Build the 2nd bubble polyline
  surface_polyline_pt[1] =
    new TriangleMeshPolyLine(bubble_vertex, Second_bubble_boundary_id);

  // Define coordinates of a point inside the bubble
  Vector<double> bubble_center(2);
  bubble_center[0] = x_center;
  bubble_center[1] = y_center;

  surface_polyline_pt[0]->set_maximum_length(0.02);
  surface_polyline_pt[1]->set_maximum_length(0.02);

  // Create closed polygon from two polylines
  Surface_polygon_pt =
    new TriangleMeshPolygon(surface_polyline_pt, bubble_center);

  Surface_polygon_pt->set_polyline_refinement_tolerance(0.05);
  Surface_polygon_pt->set_polyline_unrefinement_tolerance(0.01);
}

template<class ELEMENT>
void RelaxingBubbleProblem<ELEMENT>::generate_mesh()
{
  generate_fluid_mesh();
  generate_surface_mesh();

  // add_sub_mesh(Fluid_mesh_pt);
  // add_sub_mesh(Surface_mesh_pt);
}

template<class ELEMENT>
void RelaxingBubbleProblem<ELEMENT>::generate_fluid_mesh()
{
  // Target area for initial mesh
  double uniform_element_area = 5e-2;

  TriangleMeshClosedCurve* outer_closed_curve_pt = Outer_boundary_polygon_pt;
  Vector<TriangleMeshClosedCurve*> surface_closed_curve_pt(1);
  surface_closed_curve_pt[0] = Surface_polygon_pt;

  // Use the TriangleMeshParameters object for gathering all
  // the necessary arguments for the TriangleMesh object
  TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

  // Define the holes on the boundary
  triangle_mesh_parameters.internal_closed_curve_pt() = surface_closed_curve_pt;

  // Define the maximum element areas
  triangle_mesh_parameters.element_area() = uniform_element_area;

  // Create the mesh
  Fluid_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
    triangle_mesh_parameters, this->time_stepper_pt());

  // Set error estimator for fluid mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Fluid_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  Fluid_mesh_pt->max_permitted_error() = 2e-5;
  Fluid_mesh_pt->min_permitted_error() = 5e-6;
  Fluid_mesh_pt->max_element_size() = 5e-2;
  Fluid_mesh_pt->min_element_size() = 1e-6;
}

template<class ELEMENT>
void RelaxingBubbleProblem<ELEMENT>::generate_surface_mesh()
{
  for (unsigned i_boundary = First_bubble_boundary_id;
       i_boundary < Second_bubble_boundary_id;
       i_boundary++)
  {
    unsigned n_element = Fluid_mesh_pt->nboundary_element(i_boundary);
    for (unsigned i_element = 0; i_element < n_element; i_element++)
    {
      ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(i_boundary, i_element));

      // Find the index of the face of element e along boundary b
      int face_index =
        Fluid_mesh_pt->face_index_at_boundary(i_boundary, i_element);

      HeleShawInterfaceElement<ELEMENT>* interface_element_pt =
        new HeleShawInterfaceElement<ELEMENT>(fluid_element_pt, face_index);

      // Add it to the mesh
      Surface_mesh_pt->add_element_pt(interface_element_pt);

      // Add the appropriate boundary number
      interface_element_pt->set_boundary_number_in_bulk_mesh(i_boundary);

      interface_element_pt->q_inv_pt() = relaxing_bubble::q_inv_pt;
      interface_element_pt->st_pt() = relaxing_bubble::st_pt;
      interface_element_pt->alpha_pt() = relaxing_bubble::alpha_pt;
      interface_element_pt->upper_wall_fct_pt() =
        relaxing_bubble::upper_wall_fct;
      interface_element_pt->wall_speed_fct_pt() =
        relaxing_bubble::wall_speed_fct;
      interface_element_pt->bubble_pressure_fct_pt() =
        relaxing_bubble::bubble_pressure_fct;
    }
  }
}

#endif
