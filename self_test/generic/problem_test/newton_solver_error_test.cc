#include "generic.h"

using namespace oomph;

class SingularJacobianElement : public GeneralisedElement
{
public:
  SingularJacobianElement() : GeneralisedElement()
  {
    add_internal_data(new Data(1));
  }

  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    int local_eqn = internal_local_eqn(0, 0);
    residuals[local_eqn] = 1.0;
    // Do nothing
  }
};

class UnsolveableElement : public GeneralisedElement
{
public:
  UnsolveableElement() : GeneralisedElement()
  {
    add_internal_data(new Data(1));
  }

  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    int local_eqn = internal_local_eqn(0, 0);
    residuals[local_eqn] = 1.0;
    // Do nothing
  }

  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    int local_eqn = internal_local_eqn(0, 0);
    jacobian(local_eqn, local_eqn) = 1.0;
  }
};

void test_singular_jacobian()
{
  Mesh* mesh = new Mesh;
  mesh->add_element_pt(new SingularJacobianElement);
  Problem problem;
  problem.add_sub_mesh(mesh);
  problem.build_global_mesh();
  problem.assign_eqn_numbers();
  try
  {
    problem.newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }
}

void test_unsolveable_problem()
{
  Mesh* mesh = new Mesh;
  mesh->add_element_pt(new UnsolveableElement);
  Problem problem;
  problem.add_sub_mesh(mesh);
  problem.build_global_mesh();
  problem.assign_eqn_numbers();
  try
  {
    problem.newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }
}


int main()
{
  test_singular_jacobian();
  test_unsolveable_problem();
  return 0;
}
