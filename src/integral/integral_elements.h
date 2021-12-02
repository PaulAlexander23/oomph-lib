#ifndef INTEGRAL_ELEMENTS_HEADER
#define INTEGRAL_ELEMENTS_HEADER

#include "generic.h"

namespace oomph
{
  class IntegralEquations : public virtual FiniteElement
  {
  public:
    typedef void (*IntegrandFctPt)(const Vector<double>& x, double& f);

  private:
    IntegrandFctPt Integrand_fct_pt;

  public:
    IntegralEquations() {}

    IntegralEquations(Data* const& data_pt);

    ~IntegralEquations() {}

    void set_external_data_pt(Data* const& data_pt);

    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      bool compute_jacobian = false;
      fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, compute_jacobian);
    }

    inline virtual void get_integrand(const unsigned& ipt,
                                      const Vector<double>& x,
                                      double& f) const
    {
      if (Integrand_fct_pt == 0)
      {
        f = 0.0;
      }
      else
      {
        (*Integrand_fct_pt)(x, f);
      }
    }

    /// Access function: Pointer to integrand function
    IntegrandFctPt& integrand_fct_pt()
    {
      return Integrand_fct_pt;
    }

    /// Access function: Pointer to integrand function. Const version
    IntegrandFctPt integrand_fct_pt() const
    {
      return Integrand_fct_pt;
    }

    /// Self-test: Return 0 for OK
    unsigned self_test()
    {
      return 0;
    }

    /// Get flux: flux[i] = du/dx_i
    void get_flux(const Vector<double>& s, Vector<double>& flux) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      unsigned DIM = 2;
      // Vector for coordinates
      Vector<double> x(DIM);

      // Get x position as Vector
      interpolated_x(s, x);

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, DIM);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidx);

      // Initialise to zero
      for (unsigned j = 0; j < DIM; j++)
      {
        flux[j] = 0.0;
      }

      double f;
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < DIM; j++)
        {
          (*Integrand_fct_pt)(x, f);
          flux[j] += f * dpsidx(l, j);
        }
      }
    }

  private:
    void fill_in_generic_residual_contribution(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian,
                                               const unsigned& flag);
  };

  IntegralEquations::IntegralEquations(Data* const& data_pt)
    : GeneralisedElement()
  {
    set_external_data_pt(data_pt);
  }

  void IntegralEquations::set_external_data_pt(Data* const& data_pt)
  {
    bool use_fd_for_jacobian = true;
    this->add_external_data(data_pt, use_fd_for_jacobian);
  }


  void IntegralEquations::fill_in_generic_residual_contribution(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    const unsigned& flag)
  {
    if (this->nexternal_data() != 1)
    {
      throw OomphLibError("Should only be one external data",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    const unsigned DIM = 2;

    // Find out how many nodes there are
    const unsigned n_node = nnode();

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsidx(n_node, DIM);

    // Set the value of n_intpt
    const unsigned n_intpt = integral_pt()->nweight();

    // Integers to store the local equation and unknown numbers
    int local_eqn = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and functions
      double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      Vector<double> interpolated_x(DIM, 0.0);
      // Calculate function value
      //-----------------------------------------
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
        }
      }

      // Get source function
      //-------------------
      double f;
      get_integrand(ipt, interpolated_x, f);

      // Assemble residuals and Jacobian
      //--------------------------------

      unsigned i_data = 0;
      unsigned i_value = 0;

      // Get the local equation
      local_eqn = external_local_eqn(i_data, i_value);

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Add body force/source term here
        residuals[local_eqn] += f * psi(l) * W;
      }
    }
  }

  template<unsigned NNODE_1D>
  class QIntegralElement : public virtual QElement<2, NNODE_1D>,
                           public virtual IntegralEquations
  {
  private:
    /// Static int that holds the number of variables at
    /// nodes: always the same
    static const unsigned Initial_Nvalue;

  public:
    /// Constructor: Call constructors for QElement and
    /// Poisson equations
    QIntegralElement() : QElement<2, NNODE_1D>(), IntegralEquations() {}

    ///  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }
  };

  template<unsigned NNODE_1D>
  const unsigned QIntegralElement<NNODE_1D>::Initial_Nvalue = 0;
} // namespace oomph
#endif
