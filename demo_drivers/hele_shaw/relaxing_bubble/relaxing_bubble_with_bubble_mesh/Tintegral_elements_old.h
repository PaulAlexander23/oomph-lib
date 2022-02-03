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

  protected:
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

  private:
    /// Static unsigned that holds the (same) number of variables at every node
    static const unsigned Initial_Nvalue;
  };

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
  // Set the data for the number of Variables at each node, always 0
  //======================================================================
  template<unsigned NNODE_1D>
  const unsigned TIntegralElement<NNODE_1D>::Initial_Nvalue = 0;

  //======================================================================
  template<class TINTEGRAL_ELEMENT>
  class ProjectableTIntegralElement
    : public virtual ProjectableElement<TINTEGRAL_ELEMENT>
  {
  public:
    /// \short Specify the values associated with field fld.
    /// The information is returned in a vector of pairs which comprise
    /// the Data object and the value within it, that correspond to field fld.
    Vector<std::pair<Data*, unsigned>> data_values_of_field(const unsigned& fld)
    {
      cout << "data_values_of_field" << endl;
      // Create the vector
      Vector<std::pair<Data*, unsigned>> data_values;

      // Loop over all nodes
      unsigned nnod = this->nnode();
      for (unsigned j = 0; j < nnod; j++)
      {
        // Add the data value associated field: The node itself
        data_values.push_back(std::make_pair(this->node_pt(j), fld));
      }

      // Return the vector
      return data_values;
    }

    /// \short Number of fields to be projected: Zero
    unsigned nfields_for_projection()
    {
      cout << "nfields_for_projection" << endl;
      return 0;
    }

    /// \short Number of history values to be stored for fld-th field.
    unsigned nhistory_values_for_projection(const unsigned& fld)
    {
      cout << "nhistory_values_for_projection" << endl;
      return this->node_pt(0)->ntstorage();
    }

    ///\short Number of positional history values
    unsigned nhistory_values_for_coordinate_projection()
    {
      cout << "nhistory_values_for_coordinate_projection" << endl;
      return this->node_pt(0)->ntstorage();
    }

    /// \short Return Jacobian of mapping and shape functions of field fld
    /// at local coordinate s
    double jacobian_and_shape_of_field(const unsigned& fld,
                                       const Vector<double>& s,
                                       Shape& psi)
    {
      cout << "jacobian_and_shape_of_field" << endl;
      return 0;
    }

    /// \short Return interpolated field fld at local coordinate s, at time
    /// level t (t=0: present; t>0: history values)
    double get_field(const unsigned& t,
                     const unsigned& fld,
                     const Vector<double>& s)
    {
      cout << "get_field" << endl;
      return 0;
    }

    /// Return number of values in field fld: One per node
    unsigned nvalue_of_field(const unsigned& fld)
    {
      cout << "nvalue_of_field" << endl;
      return 0;
    }

    /// Return local equation number of value j in field fld.
    int local_equation(const unsigned& fld, const unsigned& j)
    {
      cout << "local_equation" << endl;
      return 0;
    }
  };

  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<ProjectableTIntegralElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };
} // namespace oomph
#endif
