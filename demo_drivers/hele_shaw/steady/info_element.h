#ifndef DATA_ELEMENT_HEADER
#define DATA_ELEMENT_HEADER

#include "generic.h"

namespace oomph
{
  class InfoElement : public virtual GeneralisedElement
  {
  public:
    InfoElement();

    InfoElement(Data* const& data_pt);

    ~InfoElement() {}

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution(residuals, jacobian, 1);
      GeneralisedElement::fill_in_jacobian_from_internal_by_fd(
        residuals, jacobian, false); /// Alice
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      //            std::cout << "Calling jacobian and mass matrix in flux " <<
      //            std::endl;
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }

  private:
    void fill_in_generic_residual_contribution(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian,
                                               const unsigned& flag);
  };

  InfoElement::InfoElement() : GeneralisedElement() {}

  InfoElement::InfoElement(Data* const& data_pt) : GeneralisedElement()
  {
    bool use_fd_for_jacobian = true;
    this->add_internal_data(data_pt, use_fd_for_jacobian);
  }


  void InfoElement::fill_in_generic_residual_contribution(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    const unsigned& flag)
  {
    unsigned n_data = this->ninternal_data();
    for (unsigned i_data = 0; i_data < n_data; i_data++)
    {
      unsigned n_value = this->internal_data_pt(i_data)->nvalue();
      for (unsigned i_value = 0; i_value < n_value; i_value++)
      {
        int local_eqn = this->internal_local_eqn(i_data, i_value);

        /// The contribution to the data is minus the value as it is moved from
        /// the right hand side of the equation
        residuals[local_eqn] += -this->internal_data_pt(i_data)->value(i_value);
      }
    }
  }
} // namespace oomph
#endif
