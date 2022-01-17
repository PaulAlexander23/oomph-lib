#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"

#include "filled_channel_flow_problem.h"

using namespace oomph;
using namespace std;

template<class ELEMENT>
class THeleShawChannelProblem : public HeleShawChannelProblem<ELEMENT>
{
public:
  THeleShawChannelProblem() : HeleShawChannelProblem<ELEMENT>() {}

  ~THeleShawChannelProblem() {}

  void generate_bulk_mesh();
};
template<class ELEMENT>
void THeleShawChannelProblem<ELEMENT>::generate_bulk_mesh()
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
  TriangleMeshPolygon* Outer_boundary_polygon_pt = new TriangleMeshPolygon(boundary_polyline_pt);

  // Convert to "closed curve" objects
  TriangleMeshClosedCurve* rect_closed_curve_pt = Outer_boundary_polygon_pt;

  // Generate mesh parameters for external mesh generator "Triangle"
  TriangleMeshParameters triangle_mesh_parameters(rect_closed_curve_pt);

  double maximum_default_element_area = 1e-2;
  triangle_mesh_parameters.element_area() = maximum_default_element_area;

  // Call external mesh generator
  this->Bulk_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);
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
  THeleShawChannelProblem<THeleShawElement<3>> problem;

  problem.setup();

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
