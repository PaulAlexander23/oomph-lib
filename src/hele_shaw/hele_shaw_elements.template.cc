#include "hele_shaw_elements.h"

#ifndef OOMPH_HELE_SHAW_ELEMENTS_TEMPLATE
#define OOMPH_HELE_SHAW_ELEMENTS_TEMPLATE

//======================================================================
/// QHeleShawElement elements are linear/quadrilateral/brick-shaped
/// HeleShaw elements with isoparametric interpolation for the function.
//======================================================================
// JACK - SPECIFIC ELEMENT

/// \short Output function:
///  x,y,u   or    x,y,z,u
template<unsigned NNODE_1D>
void QHeleShawElement<NNODE_1D>::output(std::ostream& outfile)
{
  HeleShawEquations::output(outfile);
}

///  \short Output function:
///   x,y,u   or    x,y,z,u at n_plot^2 plot points
template<unsigned NNODE_1D>
void QHeleShawElement<NNODE_1D>::output(std::ostream& outfile,
                                        const unsigned& n_plot)
{
  HeleShawEquations::output(outfile, n_plot);
}

/// \short C-style output function:
///  x,y,u   or    x,y,z,u
template<unsigned NNODE_1D>
void QHeleShawElement<NNODE_1D>::output(FILE* file_pt)
{
  HeleShawEquations::output(file_pt);
}

///  \short C-style output function:
///   x,y,u   or    x,y,z,u at n_plot^2 plot points
template<unsigned NNODE_1D>
void QHeleShawElement<NNODE_1D>::output(FILE* file_pt, const unsigned& n_plot)
{
  HeleShawEquations::output(file_pt, n_plot);
}

/// \short Output function for an exact solution:
///  x,y,u_exact   or    x,y,z,u_exact at n_plot^2 plot points
template<unsigned NNODE_1D>
void QHeleShawElement<NNODE_1D>::output_fct(
  std::ostream& outfile,
  const unsigned& n_plot,
  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
  HeleShawEquations::output_fct(outfile, n_plot, exact_soln_pt);
}

/// \short Output function for a time-dependent exact solution.
///  x,y,u_exact   or    x,y,z,u_exact at n_plot^2 plot points
/// (Calls the steady version)
template<unsigned NNODE_1D>
void QHeleShawElement<NNODE_1D>::output_fct(
  std::ostream& outfile,
  const unsigned& n_plot,
  const double& time,
  FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
{
  HeleShawEquations::output_fct(outfile, n_plot, time, exact_soln_pt);
}


// Inline functions:

//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned NNODE_1D>
double QHeleShawElement<NNODE_1D>::dshape_and_dtest_eulerian_hele_shaw(
  const Vector<double>& s,
  Shape& psi,
  DShape& dpsidx,
  Shape& test,
  DShape& dtestdx) const
{
  // Call the geometrical shape functions and derivatives
  const double J = this->dshape_eulerian(s, psi, dpsidx);

  // Set the test functions equal to the shape functions
  test = psi;
  dtestdx = dpsidx;

  // Return the jacobian
  return J;
}

//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned NNODE_1D>
double QHeleShawElement<NNODE_1D>::dshape_and_dtest_eulerian_at_knot_hele_shaw(
  const unsigned& ipt, Shape& psi, DShape& dpsidx, Shape& test, DShape& dtestdx)
  const
{
  // Call the geometrical shape functions and derivatives
  const double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

  // Set the pointers of the test functions
  test = psi;
  dtestdx = dpsidx;

  // Return the jacobian
  return J;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/// Define the shape functions (psi) and test functions (test) and
/// their derivatives w.r.t. global coordinates (dpsidx and dtestdx)
/// and return Jacobian of mapping (J). Additionally compute the
/// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned NNODE_1D>
double QHeleShawElement<NNODE_1D>::dshape_and_dtest_eulerian_at_knot_hele_shaw(
  const unsigned& ipt,
  Shape& psi,
  DShape& dpsidx,
  RankFourTensor<double>& d_dpsidx_dX,
  Shape& test,
  DShape& dtestdx,
  RankFourTensor<double>& d_dtestdx_dX,
  DenseMatrix<double>& djacobian_dX) const
{
  // Call the geometrical shape functions and derivatives
  const double J =
    this->dshape_eulerian_at_knot(ipt, psi, dpsidx, djacobian_dX, d_dpsidx_dX);

  // Set the pointers of the test functions
  test = psi;
  dtestdx = dpsidx;
  d_dtestdx_dX = d_dpsidx_dX;

  // Return the jacobian
  return J;
}

//=======================================================================
/// Face geometry for the QHeleShawElement elements: The spatial
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<QHeleShawElement<NNODE_1D>>
  : public virtual QElement<1, NNODE_1D>
{
public:
  /// \short Constructor: Call the constructor for the
  /// appropriate lower-dimensional QElement
  FaceGeometry() : QElement<1, NNODE_1D>() {}
};

//======================================================================
/// Set the data for the number of Variables at each node, always one
/// in every case
//======================================================================
template<unsigned NNODE_1D>
const unsigned QHeleShawElement<NNODE_1D>::Initial_Nvalue = 1;


//====================================================================
/// Force build of templates
//====================================================================
// template class QHeleShawElement<2>;
// template class QHeleShawElement<3>;
// template class QHeleShawElement<4>;

#endif
