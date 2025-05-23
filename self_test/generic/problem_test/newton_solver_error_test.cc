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

template<class ELEMENT>
Problem* create_single_element_problem()
{
  Mesh* mesh = new Mesh;
  mesh->add_element_pt(new UnsolveableElement);
  Problem* problem_pt = new Problem;
  problem_pt->add_sub_mesh(mesh);
  problem_pt->build_global_mesh();
  problem_pt->assign_eqn_numbers();
  return problem_pt;
}

void test_singular_jacobian()
{
  Problem* problem_pt =
    create_single_element_problem<SingularJacobianElement>();
  try
  {
    problem_pt->newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }
}


void test_unsolveable_problem()
{
  Problem* problem_pt = create_single_element_problem<UnsolveableElement>();
  try
  {
    problem_pt->newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }
}

void test_handling_of_multiple_exceptions()
{
  Problem* problem_pt = create_single_element_problem<UnsolveableElement>();
  try
  {
    problem_pt->newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }
  try
  {
    problem_pt->newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }
}

void test_handling_of_exceptions_with_multiple_problems()
{
  Problem* problem_pt = create_single_element_problem<UnsolveableElement>();
  Problem* problem2_pt = create_single_element_problem<UnsolveableElement>();
  try
  {
    problem_pt->newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }
  try
  {
    problem2_pt->newton_solve();
  }
  catch (OomphLibError& err)
  {
    cout << "Caught OomphLibError: " << err.what() << endl;
  }

  delete problem_pt;
  problem_pt = 0;

  try
  {
    problem2_pt->newton_solve();
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
  test_handling_of_multiple_exceptions();
  test_handling_of_exceptions_with_multiple_problems();
  return 0;
}
