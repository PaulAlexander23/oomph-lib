
#include <iostream>
#include "generic.h"
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"
#include "solid.h"
#include "constitutive.h"

using namespace std;
using namespace oomph;

int main()
{
  cout << "QTaylorHoodElements Test" << endl;


  Vector<SolidNode*> nodes;

  cout << "Axisym element" << endl;
  // AxisymmetricQTaylorHoodElement element;
  // AxisymmetricQCrouzeixRaviartElement element;
  PseudoSolidNodeUpdateElement<AxisymmetricTTaylorHoodElement,
                               TPVDElement<2, 3>>
    element;
  // AxisymmetricTCrouzeixRaviartElement element;

  // cout << element.ndof() << endl;
  cout << "add nodes" << endl;
  Vector<double> s;
  for (unsigned n = 0; n < element.nnode(); n++)
  {
    nodes.push_back(dynamic_cast<SolidNode*>(element.construct_node(n)));
    element.local_coordinate_of_node(n, s);
    for (unsigned i = 0; i < 2; i++)
    {
      nodes[n]->x(i) = s[i];
    }
  }

  cout << nodes[5]->xi(0) << endl;
  cout << nodes[5]->xi(1) << endl;
  element.node_update();
  cout << nodes[5]->xi(0) << endl;
  cout << nodes[5]->xi(1) << endl;

  Steady<1> timestepper;
  Time my_time;
  for (unsigned n = 0; n < element.nnode(); n++)
  {
    nodes[n]->time_stepper_pt() = &timestepper;
    nodes[n]->time_stepper_pt()->time_pt() = &my_time;
  }
  cout << "time: " << nodes[0]->time_stepper_pt()->time_pt()->time();

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

  element.constitutive_law_pt() = new GeneralisedHookean(new double(0.25));
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

  // DenseMatrix<double> fd_jacobian(global_number, global_number, 0.0);
  // element.fill_in_jacobian_from_internal_by_fd(residuals, fd_jacobian);
  // fd_jacobian.output(cout);

  ofstream outfile;
  outfile.open("test.dat");
  element.output(outfile);
  outfile.close();

  return (EXIT_SUCCESS);
}
