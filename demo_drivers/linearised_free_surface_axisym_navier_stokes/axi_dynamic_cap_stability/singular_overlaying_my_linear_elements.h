#ifndef SINGULAR_OVERLAYING_MY_LINEAR_ELEMENT_HEADER
#define SINGULAR_OVERLAYING_MY_LINEAR_ELEMENT_HEADER

#include "overlaying_my_linear_element.h"

namespace oomph
{
  template<class BASE_ELEMENT>
  class SingularOverlayingMyLinearElement
    : public virtual OverlayingMyLinearElement<BASE_ELEMENT>
  {
  private:
    /// Flag to indicate if the element has been augmented
    bool IsAugmented;

    /// Flag to indicate if the Jacobian is computed using finite differences
    bool IsJacobianFD;

    /// Vector of pointers to SingularNavierStokesSolutionElement objects
    Vector<
      SingularNavierStokesSolutionElement<OverlayingMyLinearElement<BASE_ELEMENT>>*>
      C_equation_elements_pt;

  public:
    SingularOverlayingMyLinearElement()
          : OverlayingMyLinearElement<BASE_ELEMENT>(),
        IsAugmented(false),
        IsJacobianFD(false)
    {
    }

    bool is_augmented() const
    {
      return IsAugmented;
    }

    // Check if the element is using finite differences for the Jacobian
    bool is_using_fd_jacobian()
    {
      return IsJacobianFD;
    }

    void augment()
    {
      IsAugmented = true;
    }

    void add_c_equation_element_pt(
      SingularNavierStokesSolutionElement<
        OverlayingMyLinearElement<BASE_ELEMENT>>* c_pt)
    {
      // Add the element
      C_equation_elements_pt.push_back(c_pt);

      // Add the additional unknown of this object as external data in the
      // Navier-Stokes element
      const bool use_fd = false;
      this->add_external_data(c_pt->internal_data_pt(0), use_fd);
    }

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      // Get the contribution from the underlying wrapped element first
      OverlayingMyLinearElement<
        BASE_ELEMENT>::fill_in_contribution_to_residuals(residuals);

      if (this->is_augmented())
      {
        this->fill_in_generic_residual_contribution_wrapped_axi_nst(
          residuals, GeneralisedElement::Dummy_matrix, 0);
      }
    }

    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Use finite differences
      if (this->is_using_fd_jacobian())
      {
        FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      }
      // Otherwise use analytic contributions
      else
      {
        // Call the base fill_in_contribution_to_jacobian function
        OverlayingMyLinearElement<BASE_ELEMENT>::fill_in_contribution_to_jacobian(
          residuals, jacobian);
        // Then call the singular Navier-Stokes element's
        // fill_in_contribution_to_jacobian function
        if (this->is_augmented())
        {
          this->fill_in_generic_residual_contribution_wrapped_axi_nst(
            residuals, jacobian, 1);
        }
      }
    }

    /// Fill in the additional contributions to the momentum equations and
    /// implement the total velocity equations
    void fill_in_generic_residual_contribution_wrapped_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
    }
  };

  template<class BASE_ELEMENT>
  class FaceGeometry<SingularOverlayingMyLinearElement<BASE_ELEMENT>>
    : public TElement<1, 3>
  {
  };
}; // namespace oomph

#endif
