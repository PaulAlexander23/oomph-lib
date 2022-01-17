#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"

#include "filled_channel_flow_problem.h"

using namespace oomph;
using namespace std;

template<class ELEMENT>
class QHeleShawChannelProblem : public HeleShawChannelProblem<ELEMENT>
{
public:
  QHeleShawChannelProblem() : HeleShawChannelProblem<ELEMENT>() {}

  ~QHeleShawChannelProblem(){};

  void generate_bulk_mesh();
};

template<class ELEMENT>
void QHeleShawChannelProblem<ELEMENT>::generate_bulk_mesh()
{
  unsigned n_x = 10;
  unsigned n_y = 10;
  double l_x = problem_parameter::domain_length;
  double l_y = problem_parameter::domain_width;

  cout << "SimpleRectangularQuadMesh" << endl;
  this->Bulk_mesh_pt =
    new SimpleRectangularQuadMesh<ELEMENT>(n_x, n_y, l_x, l_y);
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
  QHeleShawChannelProblem<QHeleShawElement<3>> problem;

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
