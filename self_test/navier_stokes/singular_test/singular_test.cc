
#include <iostream>
#include "generic.h"
#include "navier_stokes.h"
#include "solid.h"
#include "my_navier_stokes_elements_with_singularity.h"

using namespace std;
using namespace oomph;

int main()
{
  cout << "Pseudo-solid node update T Taylor Hood PVD Traction Element Test"
       << endl;

  typedef NavierStokesElementWithSingularity<
    PseudoSolidNodeUpdateElement<TTaylorHoodElement<2>, TPVDElement<2, 3>>>
    ELEMENT;

  ELEMENT* element_pt = new ELEMENT;

  Vector<Node*> nodes;
  Vector<double> s;
  for (unsigned n = 0; n < element_pt->nnode(); n++)
  {
    nodes.push_back(element_pt->construct_boundary_node(n));
    element_pt->local_coordinate_of_node(n, s);
    for (unsigned i = 0; i < 2; i++)
    {
      nodes[n]->x(i) = s[i];
    }
  }

  Data* Singular_scaling_data_pt = new Data(1);

  int face_index = 0;
  // NavierStokesTractionElement<ELEMENT> face_element(
  //  element_pt, face_index, false);
  SingularNavierStokesTractionElement<ELEMENT> face_element(
    element_pt, face_index, Singular_scaling_data_pt, false);

  unsigned long global_number = 0;
  Vector<double*> Dof_pt;

  // Loop over the nodes and call their assigment functions
  unsigned long nnod = nodes.size();

  for (unsigned long i = 0; i < nnod; i++)
  {
    nodes[i]->assign_eqn_numbers(global_number, Dof_pt);
  }

  element_pt->assign_internal_eqn_numbers(global_number, Dof_pt);
  element_pt->assign_local_eqn_numbers(true);

  cout << element_pt->local_eqn_number(3) << endl;
  cout << element_pt->u_nst(3, 0) << endl;
  cout << element_pt->momentum_local_eqn(3, 0) << endl;
  cout << element_pt->u_local_unknown(3, 0) << endl;

  ofstream outfile;
  outfile.open("test.dat");
  element_pt->output(outfile, 3);
  outfile.close();

  outfile.open("test_face.dat");
  // face_element.output(outfile);
  outfile.close();

  cout << face_element.nst_u(0, 0) << endl;
  cout << face_element.nst_momentum_local_eqn(0, 0) << endl;
  cout << face_element.nst_u_local_unknown(0, 0) << endl;
  cout << face_element.nst_continuity_local_eqn(0) << endl;
  cout << face_element.interpolated_u(s, 0) << endl;
  cout << face_element.interpolated_p(s) << endl;
  unsigned dofs = Dof_pt.size();
  cout << "dof: " << dofs << endl;

  Vector<double> residuals(dofs, 0.0);
  DenseMatrix<double> jacobian(dofs, dofs, 0.0);

  cout << jacobian.ncol() << endl;

  element_pt->get_jacobian(residuals, jacobian);

  cout << jacobian.ncol() << endl;
  jacobian.output(cout);

  return (EXIT_SUCCESS);
}