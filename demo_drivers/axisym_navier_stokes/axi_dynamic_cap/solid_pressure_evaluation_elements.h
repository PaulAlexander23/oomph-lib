#ifndef SOLID_PRESSURE_EVALUATION_ELEMENTS_HEADER
#define SOLID_PRESSURE_EVALUATION_ELEMENTS_HEADER

#include "generic.h"
#include "navier_stokes.h"

namespace oomph
{
  template<class ELEMENT>
  class SolidPressureEvaluationElement
    : public virtual PressureEvaluationElement<ELEMENT>,
      public virtual SolidFiniteElement
  {
  public:
    SolidPressureEvaluationElement(FiniteElement* const& element_pt,
                                   const int& face_index,
                                   Node* const& node_pt,
                                   const unsigned& pressure_value_index)
      : PressureEvaluationElement<ELEMENT>(
          element_pt, face_index, node_pt, pressure_value_index),
        SolidFiniteElement()
    {
    }

    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default
    /// This final over-ride is required because both SolidFiniteElements
    /// and FaceElements overload zeta_nodal
    virtual double zeta_nodal(const unsigned& n,
                              const unsigned& k,
                              const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    /// Set pointer to MacroElement -- overloads generic version
    /// and uses the MacroElement
    /// also as the default for the "undeformed" configuration.
    /// This assignment must be overwritten with
    /// set_undeformed_macro_elem_pt(...) if the deformation of
    /// the solid body is driven by a deformation of the
    /// "current" Domain/MacroElement representation of it's boundary.
    /// Can be overloaded in derived classes to perform additional
    /// tasks
    virtual void set_macro_elem_pt(MacroElement* macro_elem_pt)
    {
      Macro_elem_pt = macro_elem_pt;
      Undeformed_macro_elem_pt = macro_elem_pt;
    }

    /// Fill in contribution from Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      PressureEvaluationElement<ELEMENT>::fill_in_contribution_to_jacobian(
        residuals, jacobian);

      // Get the solid entries in the jacobian using finite differences
      fill_in_jacobian_from_solid_position_by_fd(residuals, jacobian);
    }
  };

}; // namespace oomph

#endif
