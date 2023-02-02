
#include <iostream>
#include "generic.h"

using namespace std;
using namespace oomph;

int main()
{
  Data* data_pt = new Data(2);
  unsigned long global_number = 0;
  Vector<double*> dof_pt;
  data_pt->assign_eqn_numbers(global_number, dof_pt);
  cout << data_pt->eqn_number(1) << endl;

  Data* new_data_pt = new Data(2);
  new_data_pt->assign_eqn_numbers(global_number, dof_pt);

  GeneralisedElement element;
  element.add_external_data(new_data_pt);

  element.assign_internal_eqn_numbers(global_number, dof_pt);
  element.assign_local_eqn_numbers(false);

  cout << element.local_eqn_number(0) << endl;
  cout << element.local_eqn_number(1) << endl;
  cout << element.eqn_number(0) << endl;
  cout << element.eqn_number(1) << endl;

  return 0;
}
