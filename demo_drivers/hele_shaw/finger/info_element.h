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

    void add_data_pt(Data* const& data_pt);

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      bool compute_jacobian = false;
      fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, compute_jacobian);
    }

  private:
    void fill_in_generic_residual_contribution(Vector<double>& residuals,
                                               DenseMatrix<double>& jacobian,
                                               const unsigned& flag);
  };

  InfoElement::InfoElement() : GeneralisedElement() {}

  InfoElement::InfoElement(Data* const& data_pt) : GeneralisedElement()
  {
    this->add_data_pt(data_pt);
  }

  void InfoElement::add_data_pt(Data* const& data_pt)
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
        if (local_eqn >= 0)
        {
          /// The contribution to the data is minus the value as it is moved
          /// from the right hand side of the equation
          cout << "Here!" << endl;
          residuals[local_eqn] +=
            -this->internal_data_pt(i_data)->value(i_value) + 1234567;
        }
      }
    }
  }
} // namespace oomph
#endif
