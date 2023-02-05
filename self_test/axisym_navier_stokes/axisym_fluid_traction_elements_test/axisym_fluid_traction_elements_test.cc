
#include <iostream>
#include "generic.h"
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"

using namespace std;
using namespace oomph;

int main()
{
  cout << "QTaylorHoodElements Test" << endl;

  Vector<Node*> nodes;
  cout << "Qelement" << endl;
  QElement<2, 3> qelement;

  cout << "Axisym element" << endl;
  typedef AxisymmetricTTaylorHoodElement ELEMENT;
  // AxisymmetricQTaylorHoodElement element;
  // AxisymmetricQCrouzeixRaviartElement element;
  ELEMENT element;
  // AxisymmetricTCrouzeixRaviartElement element;

  // cout << element.ndof() << endl;
  cout << "add nodes" << endl;
  Vector<double> s;
  for (unsigned n = 0; n < element.nnode(); n++)
  {
    nodes.push_back(element.construct_node(n));
    element.local_coordinate_of_node(n, s);
    for (unsigned i = 0; i < 2; i++)
    {
      nodes[n]->x(i) = s[i];
    }
  }

  Steady<1> timestepper;
  Time my_time;
  for (unsigned n = 0; n < element.nnode(); n++)
  {
    nodes[n]->time_stepper_pt() = &timestepper;
    nodes[n]->time_stepper_pt()->time_pt() = &my_time;
  }
  cout << "time: " << nodes[0]->time_stepper_pt()->time_pt()->time();


  AxisymmetricNavierStokesTractionElement<ELEMENT> face_element;
  int face_index = 0;
  face_element.build(&element, face_index);


  cout << element.ndof() << endl;

  unsigned long global_number = 0;
  Vector<double*> Dof_pt;

  // Loop over the nodes and call their assigment functions
  unsigned long nnod = nodes.size();

  for (unsigned long i = 0; i < nnod; i++)
  {
    nodes[i]->assign_eqn_numbers(global_number, Dof_pt);
  }

  element.assign_internal_eqn_numbers(global_number, Dof_pt);
  element.assign_local_eqn_numbers(true);

  cout << element.local_eqn_number(3) << endl;
  cout << element.u_nst(3, 0) << endl;
  cout << element.momentum_local_eqn(3, 0) << endl;
  cout << element.u_local_unknown(3, 0) << endl;

  element.re_pt() = new double(1.0);
  element.re_st_pt() = new double(1.0);


  cout << "Get residuals" << endl;
  Vector<double> residuals(global_number, 0.0);
  element.get_residuals(residuals);
  for (unsigned n = 0; n < residuals.size(); n++)
  {
    cout << "n: " << n << ", res: " << residuals[n] << endl;
  }

  cout << "Get jacobian" << endl;
  DenseMatrix<double> jacobian(global_number, global_number, 0.0);
  element.get_jacobian(residuals, jacobian);

  jacobian.output(cout);

  ofstream outfile;
  outfile.open("test.dat");
  element.output(outfile);
  outfile.close();

  face_element.assign_internal_eqn_numbers(global_number, Dof_pt);
  face_element.assign_local_eqn_numbers(true);

  cout << face_element.nst_u(0, 0) << endl;
  cout << face_element.nst_momentum_local_eqn(0, 0) << endl;
  cout << face_element.nst_u_local_unknown(0, 0) << endl;
  cout << face_element.nst_continuity_local_eqn(0) << endl;
  cout << face_element.interpolated_u(s, 0) << endl;
  cout << face_element.interpolated_p(s) << endl;

  return (EXIT_SUCCESS);
}
