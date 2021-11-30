#ifndef MY_CONSTRAINT_ELEMENTS_HEADER
#define MY_CONSTRAINT_ELEMENTS_HEADER

#include "generic.h"

namespace oomph
{
  class MyConstraintElement : public GeneralisedElement
  {
  private:
    /// Pointer to the desired value of the volume
    double* Prescribed_volume_pt;

    /// Pointer to the actual value of the volume
    double* Actual_volume_pt;

    /// Index of the value in traded pressure data that corresponds to the
    /// traded pressure
    unsigned Index_of_traded_pressure_value;

    /// The Data that contains the traded pressure is stored
    /// as external or internal Data for the element. What is its index
    /// in these containers?
    unsigned External_or_internal_data_index_of_traded_pressure;

  public:
    /// Constructor: Pass pointer to target volume, pointer to Data
    /// item whose value specified by index_of_traded_pressure represents
    /// the "Pressure" value that "traded" for the volume contraint.
    /// The Data is stored as external Data for this element.
    MyConstraintElement(double* prescribed_volume_pt,
                        double* actual_volume_pt,
                        Data* p_traded_data_pt,
                        const unsigned& index_of_traded_pressure)
    {
      // Store pointer to prescribed volume
      Prescribed_volume_pt = prescribed_volume_pt;

      // Store pointer to actual volume
      Actual_volume_pt = actual_volume_pt;

      bool fd_jacobian = true;
      // Add as external data and record the index
      External_or_internal_data_index_of_traded_pressure =
        add_external_data(p_traded_data_pt, fd_jacobian);

      // Record index
      Index_of_traded_pressure_value = index_of_traded_pressure;
    }

    /// Fill in the residuals for the volume constraint
    void fill_in_generic_contribution_to_residuals_my_constraint(
      Vector<double>& residuals)
    {
      // Note: This element can only be used with the associated
      // VolumeConstraintBoundingElement elements which compute the actual
      // enclosed volume; here we only add the contribution to the
      // residual; everything else, incl. the derivatives of this
      // residual w.r.t. the nodal positions of the
      // VolumeConstraintBoundingElements
      // is handled by them
      const int local_eqn = this->external_local_eqn(
        External_or_internal_data_index_of_traded_pressure,
        Index_of_traded_pressure_value);

      if (local_eqn >= 0)
      {
        cout << "my constraint element" << endl;
        cout << local_eqn << endl;
        cout << eqn_number(local_eqn) << endl;
        cout << "Volume diff: " << -*Prescribed_volume_pt + *Actual_volume_pt
             << endl;
        residuals[local_eqn] += -*Prescribed_volume_pt + *Actual_volume_pt;
      }
    }

    /// Fill in the residuals for the volume constraint
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      this->fill_in_generic_contribution_to_residuals_my_constraint(residuals);
    }

    /// Fill in the residuals and jacobian for the volume constraint
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // No contribution to jacobian; see comment in that function
      this->fill_in_generic_contribution_to_residuals_my_constraint(residuals);
    }

    /// Fill in the residuals, jacobian and mass matrix for the volume
    /// constraint.
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // No contribution to jacobian or mass matrix; see comment in that
      // function
      this->fill_in_generic_contribution_to_residuals_my_constraint(residuals);
    }
  };

} // namespace oomph

#endif
