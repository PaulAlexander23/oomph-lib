#ifndef PROJECTABLE_AXISYMMTERIC_TTAYLOR_HOOD_ELEMENTS_HEADER
#define PROJECTABLE_AXISYMMTERIC_TTAYLOR_HOOD_ELEMENTS_HEADER

#include "solid/solid_elements.h"
#include "navier_stokes.h"
#include "fluid_interface/constrained_volume_elements.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  class AxisymmetricTTaylorHoodPVDElement
    : public virtual PseudoSolidNodeUpdateElement<
        AxisymmetricTTaylorHoodElement,
        TPVDElement<2, 3>>,
      public virtual DebugJacobianSolidFiniteElement
  {
  public:
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      AxisymmetricNavierStokesEquations::fill_in_contribution_to_jacobian(
        residuals, jacobian);

      PVDEquations<2>::fill_in_contribution_to_jacobian(residuals, jacobian);

      // Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives(jacobian);

      // Fill in the external data entries by finite difference.
      // fill_in_jacobian_from_external_by_fd(residuals, jacobian, true);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      AxisymmetricNavierStokesEquations::
        fill_in_contribution_to_jacobian_and_mass_matrix(
          residuals, jacobian, mass_matrix);

      PVDEquations<2>::fill_in_contribution_to_jacobian(residuals, jacobian);

      //   Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives(jacobian);
    }
  };


  class ProjectableAxisymmetricTTaylorHoodPVDElement
    : public virtual ElementWithError<ProjectableAxisymmetricTaylorHoodElement<
        AxisymmetricTTaylorHoodPVDElement>>,
      public virtual DebugJacobianSolidFiniteElement
  {
  public:
    void pin()
    {
      for (unsigned i = 0; i < 6; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          this->node_pt(i)->pin(j);
        }
      }
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 1; j++)
        {
          this->node_pt(i)->pin(3 + j);
        }
      }
    }

    void pin_pressure()
    {
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 1; j++)
        {
          this->node_pt(i)->pin(3 + j);
        }
      }
    }
  }; // namespace >

  // Ensure that the FaceGeometry of the new element has been set up
  template<>
  class FaceGeometry<ProjectableAxisymmetricTTaylorHoodPVDElement>
    : public virtual SolidTElement<1, 3>
  {
  public:
    FaceGeometry() : SolidTElement<1, 3>() {}
  };

  // Ensure that the FaceGeometry of the new element has been set up
  template<>
  class FaceGeometry<FaceGeometry<ProjectableAxisymmetricTTaylorHoodPVDElement>>
    : public virtual SolidPointElement
  {
  public:
    FaceGeometry() : SolidPointElement() {}
  };
} // namespace oomph

#endif
