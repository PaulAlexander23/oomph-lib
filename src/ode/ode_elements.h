// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_ODE_ELEMENTS_H
#define OOMPH_ODE_ELEMENTS_H

#include "../generic/oomph_definitions.h"
#include "../generic/oomph_utilities.h"

#include "../generic/matrices.h"
#include "../generic/Vector.h"
#include "../generic/elements.h"
#include "../generic/timesteppers.h"

namespace oomph
{
  /// Element for integrating an initial value ODE
  class ODEElement : public GeneralisedElement
  {
  private:
    /// Internal or external data?
    bool data_is_internal;

    /// Functor to provide both the ode and its solution/approximation
    SolutionFunctorBase* Exact_solution_pt;

    /// Use finite differences to calculate the Jacobian?
    bool Use_fd_jacobian;

    /// Number of stored values in the single data object
    unsigned Nvalue;

  public:
    /// Default constructor: null any pointers
    ODEElement()
    {
      Exact_solution_pt = 0;
      data_is_internal = true;
      Nvalue = 0;
      Use_fd_jacobian = false;
    }

    /// Constructor: Pass time stepper and a solution function pointer, then
    /// build the element.
    ODEElement(TimeStepper* time_stepper_pt,
               SolutionFunctorBase* exact_solution_pt)
    {
      build(time_stepper_pt, exact_solution_pt);
    }

    /// Store pointers, create external data.
    void build(TimeStepper* time_stepper_pt,
               SolutionFunctorBase* exact_solution_pt)
    {
      Exact_solution_pt = exact_solution_pt;

      Vector<double> exact = this->exact_solution(0);
      Nvalue = exact.size();

      add_internal_data(new Data(time_stepper_pt, Nvalue));
      data_is_internal = true;

      Use_fd_jacobian = false;
    }

    /// Store pointers, add external data.
    void build(TimeStepper* time_stepper_pt,
               SolutionFunctorBase* exact_solution_pt,
               Data* data_pt)
    {
      Exact_solution_pt = exact_solution_pt;

      bool preserve_existing_data = true;
      data_pt->set_time_stepper(time_stepper_pt, preserve_existing_data);

      add_external_data(data_pt);
      Nvalue = external_data_pt(0)->nvalue();
      data_is_internal = false;

      Use_fd_jacobian = false;
    }


    virtual ~ODEElement() {}

    unsigned nvalue() const
    {
      return Nvalue;
    }

    /// Get residuals
    virtual void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Get pointer to one-and-only internal data object
      Data* dat_pt;
      if (data_is_internal)
      {
        dat_pt = internal_data_pt(0);
      }
      else
      {
        dat_pt = external_data_pt(0);
      }

      // Get it's values
      Vector<double> u(Nvalue, 0.0);
      dat_pt->value(u);

      // Get time stepper
      TimeStepper* time_stepper_pt = dat_pt->time_stepper_pt();

      // Get continuous time
      double t = time_stepper_pt->time();

      int local_eqn;
      Vector<double> deriv = derivative_function(t, u);
      for (unsigned j = 0, nj = deriv.size(); j < nj; j++)
      {
        if (data_is_internal)
        {
          local_eqn = this->internal_local_eqn(0, j);
        }
        else
        {
          local_eqn = this->external_local_eqn(0, j);
        }
        if (local_eqn >= 0)
        {
          // Get dudt approximation from time stepper
          double dudt = time_stepper_pt->time_derivative(1, dat_pt, j);

          // Residual is difference between the exact derivative and our
          // time stepper's derivative estimate.
          residuals[local_eqn] = deriv[j] - dudt;
        }
      }
    }

    virtual void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                  DenseMatrix<double>& jacobian)
    {
      // Get residuals
      fill_in_contribution_to_residuals(residuals);

      if (!Use_fd_jacobian)
      {
        // Get pointer to one-and-only internal data object
        Data* dat_pt;
        if (data_is_internal)
        {
          dat_pt = internal_data_pt(0);
        }
        else
        {
          dat_pt = external_data_pt(0);
        }

        // We need jacobian of residual = f - dudt so subtract diagonal
        // (dudt)/du term.
        int local_eqn;
        int local_unknown;
        const double a = dat_pt->time_stepper_pt()->weight(1, 0);
        for (unsigned i = 0; i < Nvalue; i++)
        {
          if (data_is_internal)
          {
            local_eqn = this->internal_local_eqn(0, i);
          }
          else
          {
            local_eqn = this->external_local_eqn(0, i);
          }
          if (local_eqn >= 0)
          {
            local_unknown = local_eqn;
            jacobian(local_eqn, local_unknown) -= a;
          }
        }
      }
      else
      {
        // Use FD for jacobian
        GeneralisedElement::fill_in_jacobian_from_internal_by_fd(
          residuals, jacobian, true);
      }
    }

    virtual void fill_in_contribution_to_mass_matrix(Vector<double>& residuals,
                                                     DenseMatrix<double>& mm)
    {
      fill_in_contribution_to_residuals(residuals);
      for (unsigned j = 0, nj = Nvalue; j < nj; j++)
      {
        mm(j, j) = 1;
      }
    }

    /// Exact solution
    Vector<double> exact_solution(const double& t) const
    {
#ifdef PARANOID
      if (Exact_solution_pt == 0)
      {
        throw OomphLibError("No exact solution function",
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
      Vector<double> dummy_x;
      return (*Exact_solution_pt)(t, dummy_x);
    }

    /// Exact solution
    Vector<double> derivative_function(const double& t, const Vector<double>& u)
    {
#ifdef PARANOID
      if (Exact_solution_pt == 0)
      {
        throw OomphLibError("No derivative function",
                            OOMPH_EXCEPTION_LOCATION,
                            OOMPH_CURRENT_FUNCTION);
      }
#endif
      Vector<double> dummy_x;
      return Exact_solution_pt->derivative(t, dummy_x, u);
    }
  };

} // namespace oomph

#endif
