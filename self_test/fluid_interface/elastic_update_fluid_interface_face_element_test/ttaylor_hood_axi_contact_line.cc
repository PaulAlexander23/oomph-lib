#include <iostream>
#include "generic.h"
#include "navier_stokes.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"
#include "solid.h"
#include "constitutive.h"
#include "element_with_error.h"

using namespace std;
using namespace oomph;

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

  double sigma = 1.0;
  double st = 1.0;
  double ca = 1.0;

  cout << "Bulk elements" << endl;
  const unsigned dim = 2;
  typedef ElementWithError<Hijacked<ProjectableAxisymmetricTaylorHoodElement<
    PseudoSolidNodeUpdateElement<AxisymmetricTTaylorHoodElement,
                                 TPVDElement<2, 3>>>>>
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

  element_pt->re_invfr_pt() = new double(1.0);

  // Set the direction of gravity
  Vector<double> G(2);
  G[1] = -1.0;
  element_pt->g_pt() = &G;


  element_pt->complete_setup_of_dependencies();

  cout << "node: " << element_pt->nnode() << endl;


  cout << "Face elements" << endl;
  int face_index = 0;

  //  ElasticUpdateFluidInterfaceElement<FluidInterfaceElement,
  //                                     LineDerivatives,
  //                                     ELEMENT>
  //    face_element;
  typedef ElasticUpdateFluidInterfaceElement<FluidInterfaceElement,
                                             AxisymmetricDerivatives,
                                             ELEMENT>
    FACE_ELEMENT;

  FACE_ELEMENT face_element;

  face_element.build(element_pt, face_index);

  face_element.ca_pt() = &ca;

  cout << "node: " << face_element.nnode() << endl;

  for (unsigned n = 0; n < face_element.nnode(); n++)
  {
    cout << "n: " << n << ", ";
    cout << "x: (" << face_element.node_pt(n)->position(0) << ", ";
    cout << face_element.node_pt(n)->position(1) << ")" << endl;
  }


  cout << "concrete face face element." << endl;
  // FaceGeometry<FaceGeometry<TPVDElement<2, 3>>> line_element;
  // FaceGeometry<FACE_ELEMENT> line_element;

  SolidFluidInterfaceBoundingElement<FACE_ELEMENT> line_element;
  // ElasticPointFluidInterfaceBoundingElement<
  //  ElasticLineFluidInterfaceElement<ELEMENT>>
  //  line_element;

  cout << "node: " << line_element.nnode() << endl;


  const unsigned line_index = 1;
  line_element.build(&face_element, line_index);

  for (unsigned n = 0; n < line_element.nnode(); n++)
  {
    cout << "n: " << n << ", ";
    cout << "x: (" << line_element.node_pt(n)->position(0) << ", ";
    cout << line_element.node_pt(n)->position(1) << ")" << endl;
  }

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

  unsigned ndof = element_pt->ndof();
  {
    Vector<double> residuals(ndof, 0.0);
    element_pt->get_residuals(residuals);
    for (unsigned n = 0; n < residuals.size(); n++)
    {
      cout << "n: " << n << ", res: " << residuals[n] << endl;
    }

    element_pt->describe_dofs(cout, "");
    cout << "Get jacobian" << endl;
    DenseMatrix<double> jacobian(ndof, ndof, 0.0);
    element_pt->get_jacobian(residuals, jacobian);
    jacobian.output(cout);
    cout << jacobian.nrow() << endl;
    cout << jacobian.ncol() << endl;
  }
  ndof = face_element.ndof();

  cout << "Get residuals" << endl;
  {
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
  }

  face_element.fix_lagrange_multiplier(2, 0);

  line_element.assign_internal_eqn_numbers(global_number, Dof_pt);
  line_element.assign_local_eqn_numbers(true);

  line_element.sigma_pt() = &sigma;
  line_element.st_pt() = &st;
  line_element.ca_pt() = &ca;
  line_element.set_contact_angle_fct(&contact_angle_fct);
  line_element.wall_unit_normal_fct_pt() = &wall_unit_normal_fct;

  cout << line_element.nnode() << endl;
  cout << line_element.nst_u_index(0, 0) << endl;
  cout << line_element.nst_u_index(0, 1) << endl;
  cout << line_element.nst_p_index(0) << endl;
  cout << line_element.nst_momentum_index(0, 0) << endl;
  cout << line_element.nst_momentum_index(0, 1) << endl;
  cout << line_element.nst_continuity_index(0) << endl;
  cout << line_element.nst_u_local_unknown(0, 0) << endl;
  cout << line_element.nst_u(0, 0) << endl;
  cout << line_element.fsi_lagrange_multiplier_nodal_index(0) << endl;
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

  const unsigned line_ndof = line_element.ndof();
  cout << "Get residuals" << endl;
  Vector<double> line_residuals(line_ndof, 0.0);
  line_element.get_residuals(line_residuals);
  for (unsigned n = 0; n < line_residuals.size(); n++)
  {
    cout << "n: " << n << ", res: " << line_residuals[n] << endl;
  }

  cout << "Get jacobian" << endl;
  DenseMatrix<double> line_jacobian(line_ndof, line_ndof, 0.0);
  line_element.get_jacobian(line_residuals, line_jacobian);
  line_jacobian.output(cout);
  cout << line_jacobian.nrow() << endl;
  cout << line_jacobian.ncol() << endl;


  return (EXIT_SUCCESS);
}
