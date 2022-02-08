#ifndef TINTEGRAL_ELEMENTS_HEADER
#define TINTEGRAL_ELEMENTS_HEADER

#include "generic.h"
#include "integral_elements.h"

namespace oomph
{
  template<unsigned NNODE_1D>
  class TIntegralElement : public virtual TElement<2, NNODE_1D>,
                           public virtual IntegralEquations,
                           public virtual ElementWithZ2ErrorEstimator
  {
  public:
    // Constructor: Call constructors for TElement and
    /// Integral equations
    TIntegralElement() : TElement<2, NNODE_1D>(), IntegralEquations() {}

    ///  Access function for Nvalue: # of `values' (pinned or dofs)
    /// at node n (always returns the same value at every node, 1)
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 2;
    }

    /// Get 'flux' for Z2 error recovery:  Standard flux.from Poisson equations
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      this->get_flux(s, flux);
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return TElement<2, NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return TElement<2, NNODE_1D>::vertex_node_pt(j);
    }

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_integral(const Vector<double>& s,
                                                     Shape& psi,
                                                     DShape& dpsidx,
                                                     Shape& test,
                                                     DShape& dtestdx) const;


  private:
    /// Static unsigned that holds the (same) number of variables at every node
    static const unsigned Initial_Nvalue;
  };

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double TIntegralElement<NNODE_1D>::dshape_and_dtest_eulerian_integral(
    const Vector<double>& s,
    Shape& psi,
    DShape& dpsidx,
    Shape& test,
    DShape& dtestdx) const
  {
    unsigned n_node = this->nnode();

    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    for (unsigned i = 0; i < n_node; i++)
    {
      test[i] = psi[i];
      dtestdx(i, 0) = dpsidx(i, 0);
      dtestdx(i, 1) = dpsidx(i, 1);
    }

    // Return the jacobian
    return J;
  }

  //=======================================================================
  /// Face geometry for the TIntegralElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TIntegralElement<NNODE_1D>>
    : public virtual TElement<1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<1, NNODE_1D>() {}
  };

  //======================================================================
  // Set the data for the number of Variables at each node, always 1
  //======================================================================
  template<unsigned NNODE_1D>
  const unsigned TIntegralElement<NNODE_1D>::Initial_Nvalue = 0;

} // namespace oomph
#endif
