
#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "integral.h"

#include "info_element.h"

using namespace std;
using namespace oomph;

namespace problem_parameter
{
  void integrand_fct(const Vector<double>& x, double& f)
  {
    f = (-x[0] * (x[0] - 2.0)) + (x[1] - 0.5) * (x[1] - 0.5);
  }
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
  unsigned n_x = 10;
  unsigned n_y = 10;
  double l_x = 2.0;
  double l_y = 1.0;
  this->Integral_mesh_pt =
    new SimpleRectangularQuadMesh<ELEMENT>(n_x, n_y, l_x, l_y);

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

  IntegralProblem<QIntegralElement<3>> problem;

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
