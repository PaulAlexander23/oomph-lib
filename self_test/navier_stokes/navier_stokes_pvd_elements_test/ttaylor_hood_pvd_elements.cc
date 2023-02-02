
#include <iostream>
#include "generic.h"
#include "navier_stokes.h"
#include "solid.h"

using namespace std;
using namespace oomph;

int main()
{
  cout << "TTaylorHoodElements Test" << endl;

  const unsigned dim = 2;
  Vector<Node*> nodes;
  TTaylorHoodElement<dim> element;

  cout << element.ndof() << endl;
  Vector<double> s;
  for (unsigned n = 0; n < element.nnode(); n++)
  {
    nodes.push_back(element.construct_node(n));
    element.local_coordinate_of_node(n, s);
    for (unsigned i = 0; i < dim; i++)
    {
      nodes[n]->x(i) = s[i];
    }
  }

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

  ofstream outfile;
  outfile.open("test.dat");
  element.output(outfile);
  outfile.close();

  return (EXIT_SUCCESS);
}
