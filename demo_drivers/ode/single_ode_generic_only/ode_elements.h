// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#pragma once

#include <iostream>

#include "generic.h"

namespace oomph {
/// Element for solving an ODE in t
class MyODEElement : public GeneralisedElement {
public:
  typedef void (*ODEFunctionPt)(const double &t, const double &u, double &dudt);

  /// Default constructor
  MyODEElement(TimeStepper *time_stepper_pt, ODEFunctionPt ode_fct_pt) {
    this->ODE_fct_pt = ode_fct_pt;
    this->nequation = 1;

    this->add_internal_data(new Data(time_stepper_pt, this->nequation));

    this->Use_fd_jacobian = true;
  }

  virtual ~MyODEElement() {}

  /// Get residuals
  virtual void fill_in_contribution_to_residuals(Vector<double> &residuals) {
    // Get pointer to one-and-only internal data object

    unsigned data_index = 0;
    Data *dat_pt = this->internal_data_pt(data_index);

    // Get it's values
    unsigned value_index = 0;
    double u = dat_pt->value(value_index);

    // Get time stepper
    TimeStepper *time_stepper_pt = dat_pt->time_stepper_pt();

    // Get continuous time
    double t = time_stepper_pt->time();

    double deriv = 0.0;
    this->ODE_fct_pt(t, u, deriv);
    unsigned j = 0;
    // Get dudt approximation from time stepper
    double dudt = time_stepper_pt->time_derivative(1, dat_pt, j);

    // Residual is difference between the exact derivative and our
    // time stepper's derivative estimate.
    residuals[j] = deriv - dudt;
  }

  virtual void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                                DenseMatrix<double> &jacobian) {
    // Get residuals
    fill_in_contribution_to_residuals(residuals);

    if (this->Use_fd_jacobian) {
      // Use FD for jacobian
      GeneralisedElement::fill_in_jacobian_from_internal_by_fd(
          residuals, jacobian, this->Use_fd_jacobian);
    } else {
      throw OomphLibError(
          "Non finite difference jacobian use is not implemented",
          OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
    }
  }

  virtual void fill_in_contribution_to_mass_matrix(Vector<double> &residuals,
                                                   DenseMatrix<double> &mm) {
    fill_in_contribution_to_residuals(residuals);

    for (unsigned j = 0, nj = this->nequation; j < nj; j++) {
      mm(j, j) = 1;
    }
  }

  void output(std::ofstream &output_stream) {
    unsigned data_index = 0;
    double t = this->internal_data_pt(data_index)->time_stepper_pt()->time();
    unsigned value_index = 0;
    double x = this->internal_data_pt(data_index)->value(value_index);
    output_stream << "t: " << t << ", x: " << x << endl;
  }

private:
  ODEFunctionPt ODE_fct_pt;
  unsigned nequation;
  bool Use_fd_jacobian;
};
} // namespace oomph
