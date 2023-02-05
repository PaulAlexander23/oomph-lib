
#include <iostream>
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "solid.h"
#include "constitutive.h"

using namespace std;
using namespace oomph;

int main()
{
  cout << "Q Taylor Hood Face Element Test" << endl;

  const unsigned dim = 2;
  typedef PseudoSolidNodeUpdateElement<TTaylorHoodElement<2>, TPVDElement<2, 3>>
    ELEMENT;
  ELEMENT* element_pt = new ELEMENT;

  Vector<SolidNode*> nodes;
  Vector<double> s;
  for (unsigned n = 0; n < element_pt->nnode(); n++)
  {
    nodes.push_back(
      dynamic_cast<SolidNode*>(element_pt->construct_boundary_node(n)));
    element_pt->local_coordinate_of_node(n, s);
    for (unsigned i = 0; i < dim; i++)
    {
      nodes[n]->x(i) = s[i];
    }
  }
  element_pt->node_update();

  element_pt->constitutive_law_pt() = new GeneralisedHookean(new double(0.25));

  // Set the Reynolds number
  element_pt->re_pt() = new double(1.0);

  // Set the Reynolds Strouhal number
  element_pt->re_st_pt() = new double(1.0);

  element_pt->complete_setup_of_dependencies();

  cout << "here" << endl;
  cout << element_pt->nnode() << endl;


  cout << "here" << endl;
  int face_index = 0;

  //  ElasticUpdateFluidInterfaceElement<FluidInterfaceElement,
  //                                     LineDerivatives,
  //                                     ELEMENT>
  //    face_element;

  ElasticLineFluidInterfaceElement<ELEMENT> face_element;

  face_element.build(element_pt, face_index);

  face_element.ca_pt() = new double(1.0);

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

  cout << global_number << endl;
  cout << element_pt->local_eqn_number(3) << endl;
  cout << element_pt->u_nst(3, 0) << endl;
  cout << element_pt->momentum_local_eqn(3, 0) << endl;
  cout << element_pt->u_local_unknown(3, 0) << endl;

  ofstream outfile;
  outfile.open("test.dat");
  element_pt->output(outfile);
  outfile.close();

  face_element.assign_internal_eqn_numbers(global_number, Dof_pt);
  face_element.assign_local_eqn_numbers(true);

  outfile.open("test_face.dat");
  face_element.output(outfile, 3);
  outfile.close();

  cout << face_element.nst_u(0, 0) << endl;
  cout << face_element.nst_momentum_local_eqn(0, 0) << endl;
  cout << face_element.nst_u_local_unknown(0, 0) << endl;
  cout << face_element.nst_continuity_local_eqn(0) << endl;
  cout << face_element.interpolated_u(s, 0) << endl;
  cout << face_element.interpolated_p(s) << endl;
  cout << face_element.sigma(s) << endl;

  const unsigned ndof = face_element.ndof();

  cout << "Get residuals" << endl;
  Vector<double> residuals(ndof, 0.0);
  face_element.get_residuals(residuals);
  for (unsigned n = 0; n < residuals.size(); n++)
  {
    cout << "n: " << n << ", res: " << residuals[n] << endl;
  }

  cout << "Get jacobian" << endl;
  DenseMatrix<double> jacobian(ndof, ndof, 0.0);
  face_element.get_jacobian(residuals, jacobian);
  jacobian.output(cout);
  cout << jacobian.nrow() << endl;
  cout << jacobian.ncol() << endl;

  return (EXIT_SUCCESS);
}
