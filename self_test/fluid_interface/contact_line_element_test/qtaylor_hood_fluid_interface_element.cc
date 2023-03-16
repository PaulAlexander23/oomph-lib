
#include <iostream>
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"

using namespace std;
using namespace oomph;


template<class ELEMENT>
class ConcreteFluidInterfaceElement : public virtual FluidInterfaceElement,
                                      public virtual FaceGeometry<ELEMENT>
{
public:
  double compute_surface_derivatives(const Shape& psi,
                                     const DShape& dpsids,
                                     const DenseMatrix<double>& interpolated_t,
                                     const Vector<double>& interpolated_x,
                                     DShape& dpsidS,
                                     DShape& dpsidS_div)
  {
    return 1.0;
  }

  void output(FILE*) {}
  void output(FILE*, const unsigned& n_plot) {}
  void output(std::ostream& outfile) {}
  void output(std::ostream& outfile, const unsigned int& n_plot) {}
};

template<class ELEMENT>
class FaceGeometry<ConcreteFluidInterfaceElement<ELEMENT>>
  : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
{
public:
  FaceGeometry()
  {
    cout << "here" << endl;
  }
};

template<class ELEMENT>
class ConcreteContactLineElement : public virtual ContactLineElement,
                                   public virtual FaceGeometry<ELEMENT>
{
public:
  ConcreteContactLineElement<ELEMENT>()
    : ContactLineElement(), FaceGeometry<ELEMENT>()
  {
  }

  virtual const unsigned fsi_lagrange_multiplier_nodal_index(const unsigned& n)
  {
    return 0;
  }
};

double contact_angle_fct(const Vector<double>& x, const Vector<double>& local_u)
{
  return 1.0;
}

void wall_unit_normal_fct(const Vector<double>& x, Vector<double>& unit_normal)
{
  unit_normal.resize(x.size());
  unit_normal[0] = 1.0;
}

int main()
{
  cout << "Q Taylor Hood Face Element Test" << endl;

  const unsigned dim = 2;
  typedef QTaylorHoodElement<dim> ELEMENT;
  Vector<Node*> nodes;
  ELEMENT* element_pt = new ELEMENT;
  element_pt->set_n_node(9);

  cout << element_pt->ndof() << endl;
  Vector<double> s;
  for (unsigned n = 0; n < element_pt->nnode(); n++)
  {
    // nodes.push_back(element_pt->construct_node(n));
    nodes.push_back(element_pt->construct_boundary_node(n));
    element_pt->local_coordinate_of_node(n, s);
    for (unsigned i = 0; i < dim; i++)
    {
      nodes[n]->x(i) = s[i];
    }
  }

  cout << "here" << endl;
  cout << element_pt->nnode() << endl;


  cout << "here" << endl;
  int face_index = -1;
  ConcreteFluidInterfaceElement<ELEMENT> face_element;
  element_pt->build_face_element(face_index, &face_element);

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

  cout << "Bulk element type: " << endl;
  cout << element_pt->element_geometry() << endl;

  ofstream outfile;
  outfile.open("test.dat");
  element_pt->output(outfile);
  outfile.close();

  face_element.assign_internal_eqn_numbers(global_number, Dof_pt);
  face_element.assign_local_eqn_numbers(true);

  outfile.open("test_face.dat");
  // face_element.output(outfile);
  outfile.close();

  cout << face_element.nst_u(0, 0) << endl;
  cout << face_element.nst_momentum_local_eqn(0, 0) << endl;
  cout << face_element.nst_u_local_unknown(0, 0) << endl;
  cout << face_element.nst_continuity_local_eqn(0) << endl;
  cout << face_element.interpolated_u(s, 0) << endl;
  cout << face_element.interpolated_p(s) << endl;
  cout << face_element.sigma(s) << endl;

  cout << "Face element type: " << endl;
  cout << face_element.element_geometry() << endl;


  cout << "Face geometry q2 element." << endl;
  FaceGeometry<QElement<2, 2>> fg2;
  cout << fg2.element_geometry() << endl;
  cout << fg2.nnode_1d() << endl;

  cout << "Face geometry q1 element." << endl;
  FaceGeometry<QElement<1, 2>> fg4;
  cout << fg4.nnode_1d() << endl;

  // FaceGeometry<FaceGeometry<QElement<2, 2>>> fg3;
  // cout << fg3.element_geometry() << endl;
  // cout << fg3.nnode_1d() << endl;

  cout << "Face geometry face face element." << endl;
  FaceGeometry<ConcreteFluidInterfaceElement<ELEMENT>> fg;
  cout << fg.nnode_1d() << endl;

  cout << "concrete face face element." << endl;
  ConcreteContactLineElement<ConcreteFluidInterfaceElement<ELEMENT>>
    line_element;
  cout << line_element.nnode_1d() << endl;
  const unsigned line_index = 1;
  line_element.build(&face_element, line_index);

  double sigma = 1.0;
  line_element.sigma_pt() = &sigma;
  double st = 1.0;
  line_element.st_pt() = &st;
  double ca = 1.0;
  line_element.ca_pt() = &ca;
  line_element.set_contact_angle_fct(&contact_angle_fct);
  line_element.wall_unit_normal_fct_pt() = &wall_unit_normal_fct;

  line_element.assign_internal_eqn_numbers(global_number, Dof_pt);
  line_element.assign_local_eqn_numbers(true);

  cout << line_element.nnode() << endl;
  cout << line_element.nst_u_index(0, 0) << endl;
  cout << line_element.nst_u_index(0, 1) << endl;
  cout << line_element.nst_p_index(0) << endl;
  cout << line_element.nst_momentum_index(0, 0) << endl;
  cout << line_element.nst_momentum_index(0, 1) << endl;
  cout << line_element.nst_continuity_index(0) << endl;
  cout << line_element.nst_u_local_unknown(0, 0) << endl;
  cout << line_element.nst_u(0, 0) << endl;
  cout << line_element.sigma() << endl;
  cout << line_element.st() << endl;
  cout << line_element.ca() << endl;
  Vector<double> position(2, 0.0);
  Vector<double> velocity(2, 0.0);
  cout << line_element.get_contact_angle(position, velocity) << endl;
  Vector<double> wall_normal;
  line_element.wall_unit_normal(position, wall_normal);
  cout << wall_normal[0] << endl;
  cout << wall_normal[1] << endl;

  const unsigned ndof = line_element.ndof();
  cout << "Get residuals" << endl;
  Vector<double> residuals(ndof, 0.0);
  line_element.get_residuals(residuals);
  for (unsigned n = 0; n < residuals.size(); n++)
  {
    cout << "n: " << n << ", res: " << residuals[n] << endl;
  }

  cout << "Get jacobian" << endl;
  DenseMatrix<double> jacobian(ndof, ndof, 0.0);
  line_element.get_jacobian(residuals, jacobian);
  jacobian.output(cout);
  cout << jacobian.nrow() << endl;
  cout << jacobian.ncol() << endl;

  return (EXIT_SUCCESS);
}
