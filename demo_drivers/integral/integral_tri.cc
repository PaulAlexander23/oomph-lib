
#include <iostream>

#include "generic.h"
#include "meshes.h"
//#include "integral.h"

#include "integral_elements.h"
#include "Tintegral_elements.h"
#include "info_element.h"

using namespace std;
using namespace oomph;

namespace problem_parameter
{
  void integrand_fct(const Vector<double>& x, double& f)
  {
    f = (-x[0] * (x[0] - 2.0)) + (x[1] - 0.5) * (x[1] - 0.5);
  }

  enum
  {
    LOWER_BOUNDARY,
    OUTLET_BOUNDARY,
    UPPER_BOUNDARY,
    INLET_BOUNDARY
  };
}; // namespace problem_parameter

template<class ELEMENT>
class IntegralProblem : public Problem
{
public:
  /// Constructor
  IntegralProblem();

  /// Destructor
  ~IntegralProblem(){};

  /// Document the solution
  void doc_solution(DocInfo& doc_info);

private:
  /// Generate mesh
  void generate_mesh();

  /// Upcast elements and finalise setup
  void upcast_and_finalise_elements();

  /// Pointer to the outer boundary polyline
  TriangleMeshPolygon* Rect_boundary_polyline_pt;

  /// Pointer to the "bulk" integral mesh
  Mesh* Integral_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;
};

template<class ELEMENT>
IntegralProblem<ELEMENT>::IntegralProblem()
{
  cout << "Generate mesh" << endl;
  this->generate_mesh();

  cout << "Upcast and finalise elements" << endl;
  this->upcast_and_finalise_elements();

  // Setup equation numbering scheme
  cout << "Number of equations: " << endl;
  cout << this->assign_eqn_numbers() << endl;
  cout << "Number of unknowns: " << endl;
  cout << this->ndof() << endl;
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::generate_mesh()
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

  // Convert to "closed curve" objects
  TriangleMeshClosedCurve* rect_closed_curve_pt =
    this->Rect_boundary_polyline_pt;

  // Generate mesh parameters for external mesh generator "Triangle"
  TriangleMeshParameters triangle_mesh_parameters(rect_closed_curve_pt);

  double maximum_default_element_area = 1e-2;
  triangle_mesh_parameters.element_area() = maximum_default_element_area;

  // Call external mesh generator
  this->Integral_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

  this->Info_mesh_pt = new Mesh;
  unsigned n_values = 1;
  Data* integral_data_pt = new Data(n_values);
  this->Info_mesh_pt->add_element_pt(new InfoElement(integral_data_pt));

  this->add_sub_mesh(this->Integral_mesh_pt);
  this->add_sub_mesh(this->Info_mesh_pt);

  this->build_global_mesh();
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::upcast_and_finalise_elements()
{
  /// Bulk mesh
  // Find number of elements in mesh
  unsigned n_element = this->Integral_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(this->Integral_mesh_pt->element_pt(i));

    el_pt->integrand_fct_pt() = problem_parameter::integrand_fct;
    el_pt->set_external_data_pt(
      this->Info_mesh_pt->element_pt(0)->internal_data_pt(0));
  }
}

template<class ELEMENT>
void IntegralProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  string data_directory = doc_info.directory();
  string filename = data_directory + "integral.dat";
  ofstream output_stream;
  output_stream.open(filename.c_str());
  double my_integral =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value(0);
  output_stream << "Integral: " << my_integral << endl;
  output_stream.close();
}

int main(int argc, char* argv[])
{
  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Create a DocInfo object for output processing
  DocInfo doc_info;
  doc_info.set_directory("RESLT/");

  IntegralProblem<TIntegralElement<3>> problem;

  /// Run problem self test
  if (problem.self_test())
  {
    throw OomphLibError(
      "Self test failed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  DoubleVector residuals;
  DenseDoubleMatrix jacobian;
  DoubleVector residualsFD;
  DenseDoubleMatrix jacobianFD;

  problem.get_jacobian(residuals, jacobian);
  problem.get_fd_jacobian(residualsFD, jacobianFD);

  for (unsigned j = 0; j < 675; j++)
  {
    printf("j: %3u, ", j);
    for (unsigned i = 674 - 3; i <= 674; i++)
    {
      printf("act: %8.5f, exp: %8.5f,",
             jacobian(i, j),
             jacobianFD(i, j));
    }
    printf("\n");
  }


  /// Call problem solve
  problem.newton_solve();

  /// Document solution
  problem.doc_solution(doc_info);

  return 0;
}
