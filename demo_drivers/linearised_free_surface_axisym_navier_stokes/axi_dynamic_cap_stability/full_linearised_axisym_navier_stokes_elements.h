// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Header file for linearised axisymmetric Navier-Stokes elements

#ifndef OOMPH_FULL_LINEARISED_AXISYM_NAVIER_STOKES_ELEMENTS_HEADER
#define OOMPH_FULL_LINEARISED_AXISYM_NAVIER_STOKES_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib includes
#include "generic.h"

namespace oomph
{
  //=======================================================================
  /// A class for elements that solve the linearised version of the
  /// unsteady Navier--Stokes equations in cylindrical polar coordinates,
  /// where we have Fourier-decomposed in the azimuthal direction so that
  /// the theta-dependance is replaced by an azimuthal mode number.
  //=======================================================================
  class FullLinearisedAxisymmetricNavierStokesEquations
    : public virtual FiniteElement
  {
  private:
    /// Static "magic" number that indicates that the pressure is not
    /// stored at a node
    static int Pressure_not_stored_at_node;

    /// Static default value for the physical constants
    /// (all initialised to zero)
    static double Default_Physical_Constant_Value;

    /// Static default value for the azimuthal mode number (zero)
    static int Default_Azimuthal_Mode_Number_Value;

    /// Static default value for the physical ratios (all initialised to one)
    static double Default_Physical_Ratio_Value;

    /// Static default value for the gravity vector (zero)
    static Vector<double> Default_Gravity_Vector;

  protected:
    // Physical constants
    // ------------------

    /// Pointer to the viscosity ratio (relative to the
    /// viscosity used in the definition of the Reynolds number)
    double* Viscosity_Ratio_pt;

    /// Pointer to the density ratio (relative to the
    /// density used in the definition of the Reynolds number)
    double* Density_Ratio_pt;

    /// Pointer to global Reynolds number
    double* Re_pt;

    /// Pointer to global Reynolds number x Strouhal number (=Womersley)
    double* ReSt_pt;

    /// Pointer to global Reynolds number x inverse Froude number
    /// (= Bond number / Capillary number)
    double* ReInvFr_pt;

    /// Pointer to global Reynolds number x inverse Rossby number
    /// (used when in a rotating frame)
    double* ReInvRo_pt;

    /// Pointer to global gravity Vector
    Vector<double>* G_pt;

    /// Pointer to azimuthal mode number k in e^ik(theta) decomposition
    int* Azimuthal_Mode_Number_pt;

    /// Pointer to base flow solution (velocity components) function
    void (*Base_flow_u_fct_pt)(const double& time,
                               const Vector<double>& x,
                               Vector<double>& result);

    /// Pointer to derivatives of base flow solution velocity
    /// components w.r.t. global coordinates (r and z) function
    void (*Base_flow_dudx_fct_pt)(const double& time,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result);

    /// Pointer to derivatives of base flow solution velocity
    /// components w.r.t. time function
    void (*Base_flow_dudt_fct_pt)(const double& time,
                                  const Vector<double>& x,
                                  Vector<double>& result);

    /// Pointer to base flow solution (pressure) function
    void (*Base_flow_p_fct_pt)(const double& time,
                               const Vector<double>& x,
                               double& result);

    /// Pointer to derivs w.r.t. nodal coords X_{pq} of spatial
    /// derivatives of base flow solution velocities function
    void (*Base_flow_d_dudx_dX_fct_pt)(const double& time,
                                       const Vector<double>& x,
                                       RankFourTensor<double>& result);

    /// Pointer to (base flow) body force function
    void (*Body_force_fct_pt)(const double& time,
                              const Vector<double>& x,
                              Vector<double>& result);

    /// Pointer to (base flow) volumetric source function
    double (*Source_fct_pt)(const double& time, const Vector<double>& x);

    /// Pointer to derivatives of base flow solution velocity
    /// components w.r.t. local coordinates (s_1 and s_2) function
    void (*Base_flow_duds_fct_pt)(const double& time,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result);

    /// Boolean flag to indicate if ALE formulation is disabled when
    /// the time-derivatives are computed. Only set to true if you're sure
    /// that the mesh is stationary.
    bool ALE_is_disabled;

    /// Access function for the local equation number
    /// information for the i-th component of the pressure.
    /// p_local_eqn[n,i] = local equation number or < 0 if pinned.
    virtual int p_local_eqn(const unsigned& n, const unsigned& i) = 0;

    /// Compute the shape functions and their derivatives
    /// w.r.t. global coordinates at local coordinate s.
    /// Return Jacobian of mapping between local and global coordinates.
    virtual double dshape_and_dtest_eulerian_lin_axi_nst(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const = 0;

    /// Compute the shape functions and their derivatives
    /// w.r.t. global coordinates at the ipt-th integration point.
    /// Return Jacobian of mapping between local and global coordinates.
    virtual double dshape_and_dtest_eulerian_at_knot_lin_axi_nst(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const = 0;

    /// Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return Jacobian of mapping (J). Also compute
    /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
    // virtual double
    // dshape_and_dtest_eulerian_and_dnodal_coordinates_at_knot_lin_axi_nst(
    //  const unsigned& ipt,
    //  Shape& psi,
    //  DShape& dpsidx,
    //  RankFourTensor<double>& d_dpsidx_dX,
    //  Shape& test,
    //  DShape& dtestdx,
    //  RankFourTensor<double>& d_dtestdx_dX,
    //  DenseMatrix<double>& djacobian_dX) const = 0;

    /// Compute the pressure shape functions at local coordinate s
    virtual void pshape_lin_axi_nst(const Vector<double>& s,
                                    Shape& psi) const = 0;

    /// Compute the pressure shape and test functions at local coordinate s
    virtual void pshape_lin_axi_nst(const Vector<double>& s,
                                    Shape& psi,
                                    Shape& test) const = 0;

  public:
    /// Calculate the velocity components of the base flow solution
    /// at a given time and Eulerian position
    virtual void get_base_flow_u(const double& time,
                                 const unsigned& ipt,
                                 const Vector<double>& x,
                                 Vector<double>& result) const
    {
      // If the function pointer is zero return zero
      if (Base_flow_u_fct_pt == 0)
      {
        // Loop over velocity components and set base flow solution to zero
        for (unsigned i = 0; i < 3; i++)
        {
          result[i] = 0.0;
        }
      }
      // Otherwise call the function
      else
      {
        (*Base_flow_u_fct_pt)(time, x, result);
      }
    }

    virtual void get_base_flow_u(const Vector<double>& s,
                                 Vector<double>& result) const = 0;

  protected:
    /// Calculate the derivatives of the velocity components of the
    /// base flow solution w.r.t. global coordinates (r and z) at a given
    /// time and Eulerian position
    virtual void get_base_flow_dudx(const double& time,
                                    const unsigned& ipt,
                                    const Vector<double>& x,
                                    DenseMatrix<double>& result) const
    {
      // If the function pointer is zero return zero
      if (Base_flow_dudx_fct_pt == 0)
      {
        // Loop over velocity components
        for (unsigned i = 0; i < 3; i++)
        {
          // Loop over coordinate directions and set to zero
          for (unsigned j = 0; j < 2; j++)
          {
            result(i, j) = 0.0;
          }
        }
      }
      // Otherwise call the function
      else
      {
        (*Base_flow_dudx_fct_pt)(time, x, result);
      }
    }

    /// Calculate the derivative of the velocity components of the
    /// base flow solution w.r.t. time at a given time and Eulerian position
    virtual void get_base_flow_dudt(const double& time,
                                    const unsigned& ipt,
                                    const Vector<double>& x,
                                    Vector<double>& result) const
    {
      // If the function pointer is zero return zero
      if (Base_flow_dudt_fct_pt == 0)
      {
        // Loop over velocity components and set to zero
        for (unsigned i = 0; i < 3; i++)
        {
          result[i] = 0.0;
        }
      }
      // Otherwise call the function
      else
      {
        (*Base_flow_dudt_fct_pt)(time, x, result);
      }
    }

    /// Calculate the pressure in the base flow solution
    /// at a given time and Eulerian position
    virtual void get_base_flow_p(const double& time,
                                 const unsigned& ipt,
                                 const Vector<double>& x,
                                 double& result) const
    {
      // If the function pointer is zero return zero
      if (Base_flow_p_fct_pt == 0)
      {
        result = 0.0;
      }

      // Otherwise call the function
      else
      {
        (*Base_flow_p_fct_pt)(time, x, result);
      }
    }

    /// Calculate the derivatives w.r.t. nodal coordinates X_{pq} of
    /// the spatial derivatives of the velocity components of the base flow
    /// solution at a given time and Eulerian position
    virtual void get_base_flow_d_dudx_dX(const double& time,
                                         const unsigned& ipt,
                                         const Vector<double>& x,
                                         RankFourTensor<double>& result) const
    {
      // If the function pointer is zero return zero
      if (Base_flow_d_dudx_dX_fct_pt == 0)
      {
        // Determine number of nodes in this element
        const unsigned n_node = nnode();

        // Loop over the element's nodes
        for (unsigned q = 0; q < n_node; q++)
        {
          // Loop over the two coordinate directions
          for (unsigned p = 0; p < 2; p++)
          {
            // Loop over velocity components
            for (unsigned i = 0; i < 3; i++)
            {
              // Loop over coordinate directions and set to zero
              for (unsigned k = 0; k < 2; k++)
              {
                result(p, q, i, k) = 0.0;
              }
            }
          }
        }
      }
      // Otherwise call the function
      else
      {
        (*Base_flow_d_dudx_dX_fct_pt)(time, x, result);
      }
    }

    /// Calculate the body force fct of the base flow at a given
    /// time and Eulerian position
    void get_body_force_base_flow(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  Vector<double>& result)
    {
      // If the function pointer is zero return zero
      if (Body_force_fct_pt == 0)
      {
        // Loop over dimensions and set body forces to zero
        for (unsigned i = 0; i < 3; i++)
        {
          result[i] = 0.0;
        }
      }
      // Otherwise call the function
      else
      {
        (*Body_force_fct_pt)(time, x, result);
      }
    }

    /// Calculate the source fct of the base flow at given time
    /// and Eulerian position
    double get_source_base_flow(const double& time,
                                const unsigned& ipt,
                                const Vector<double>& x)
    {
      // If the function pointer is zero return zero
      if (Source_fct_pt == 0)
      {
        return 0;
      }

      // Otherwise call the function
      else
      {
        return (*Source_fct_pt)(time, x);
      }
    }

    /// Calculate the gradient of the body force of the base flow
    /// at a given time and Eulerian position
    void get_body_force_gradient_base_flow(const double& time,
                                           const unsigned& ipt,
                                           const Vector<double>& x,
                                           DenseMatrix<double>& result)
    {
      // Reference value
      Vector<double> body_force(3, 0.0);
      get_body_force_base_flow(time, ipt, x, body_force);

      // FD it
      const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;
      Vector<double> body_force_pls(3, 0.0);
      Vector<double> x_pls(x);
      for (unsigned i = 0; i < 2; i++)
      {
        x_pls[i] += eps_fd;
        get_body_force_base_flow(time, ipt, x_pls, body_force_pls);
        for (unsigned j = 0; j < 3; j++)
        {
          result(j, i) = (body_force_pls[j] - body_force[j]) / eps_fd;
        }
        x_pls[i] = x[i];
      }
    }

    /// Calculate the gradient of the source function of the base flow
    /// at a given time and Eulerian position
    void get_source_gradient_base_flow(const double& time,
                                       const unsigned& ipt,
                                       const Vector<double>& x,
                                       Vector<double>& result)
    {
      // Reference value
      const double source = get_source_base_flow(time, ipt, x);

      // FD it
      const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;
      double source_pls = 0.0;
      Vector<double> x_pls(x);
      for (unsigned i = 0; i < 2; i++)
      {
        x_pls[i] += eps_fd;
        source_pls = get_source_base_flow(time, ipt, x_pls);
        result[i] = (source_pls - source) / eps_fd;
        x_pls[i] = x[i];
      }
    }

    /// Calculate the derivatives of the velocity components of the
    /// base flow solution w.r.t. local coordinates (s_1 and s_2) at a given
    /// time and Eulerian position
    virtual void get_base_flow_duds(const double& time,
                                    const unsigned& ipt,
                                    const Vector<double>& x,
                                    DenseMatrix<double>& result) const
    {
      // If the function pointer is zero return zero
      if (Base_flow_duds_fct_pt == 0)
      {
        // Loop over velocity components
        for (unsigned i = 0; i < 3; i++)
        {
          // Loop over coordinate directions and set to zero
          for (unsigned j = 0; j < 2; j++)
          {
            result(i, j) = 0.0;
          }
        }
      }
      // Otherwise call the function
      else
      {
        (*Base_flow_duds_fct_pt)(time, x, result);
      }
    }

    virtual void get_base_flow_duds(const Vector<double>& s,
                                    DenseMatrix<double>& result) const = 0;

    virtual void get_base_flow_dudx(const Vector<double>& s,
                                    DenseMatrix<double>& result) const = 0;

    virtual void get_base_flow_dudt(const Vector<double>& s,
                                    Vector<double>& result) const = 0;

    virtual void get_base_flow_p(const Vector<double>& s,
                                 double& result) const = 0;

    /// Compute the residuals for the Navier-Stokes equations;
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    virtual void fill_in_generic_residual_contribution_lin_axi_nst(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag);

  public:
    /// Constructor: NULL the base flow solution and the
    /// derivatives of the base flow function, as well as the base flow
    /// body force and volumetric source function pointers
    FullLinearisedAxisymmetricNavierStokesEquations()
      : Base_flow_u_fct_pt(0),
        Base_flow_dudx_fct_pt(0),
        Base_flow_dudt_fct_pt(0),
        Base_flow_p_fct_pt(0),
        Base_flow_d_dudx_dX_fct_pt(0),
        Body_force_fct_pt(0),
        Source_fct_pt(0),
        Base_flow_duds_fct_pt(0),
        ALE_is_disabled(false)
    {
      // Set all the physical parameter pointers to the default value of zero
      Re_pt = &Default_Physical_Constant_Value;
      ReSt_pt = &Default_Physical_Constant_Value;
      ReInvFr_pt = &Default_Physical_Constant_Value;
      ReInvRo_pt = &Default_Physical_Constant_Value;
      G_pt = &Default_Gravity_Vector;

      // Set the azimuthal mode number to default (zero)
      Azimuthal_Mode_Number_pt = &Default_Azimuthal_Mode_Number_Value;

      // Set the physical ratios to the default value of one
      Viscosity_Ratio_pt = &Default_Physical_Ratio_Value;
      Density_Ratio_pt = &Default_Physical_Ratio_Value;
    }

    /// Vector to decide whether the stress-divergence form is used or not.
    //  N.B. This needs to be public so that the intel compiler gets things
    // correct. Somehow the access function messes things up when going to
    // refineable navier--stokes
    static Vector<double> Gamma;

    // Access functions for the physical constants
    // -------------------------------------------

    /// Reynolds number
    const double& re() const
    {
      return *Re_pt;
    }

    /// Pointer to Reynolds number
    double*& re_pt()
    {
      return Re_pt;
    }

    /// Product of Reynolds and Strouhal number (=Womersley number)
    const double& re_st() const
    {
      return *ReSt_pt;
    }

    /// Pointer to product of Reynolds and Strouhal number (=Womersley number)
    double*& re_st_pt()
    {
      return ReSt_pt;
    }

    /// Global inverse Froude number
    const double& re_invfr() const
    {
      return *ReInvFr_pt;
    }

    /// Pointer to global inverse Froude number
    double*& re_invfr_pt()
    {
      return ReInvFr_pt;
    }

    /// Global Reynolds number multiplied by inverse Rossby number
    const double& re_invro() const
    {
      return *ReInvRo_pt;
    }

    /// Pointer to global Reynolds number multiplied by inverse Rossby number
    double*& re_invro_pt()
    {
      return ReInvRo_pt;
    }

    /// Vector of gravitational components
    const Vector<double>& g() const
    {
      return *G_pt;
    }

    /// Pointer to Vector of gravitational components
    Vector<double>*& g_pt()
    {
      return G_pt;
    }

    /// Azimuthal mode number k in e^ik(theta) decomposition
    const int& azimuthal_mode_number() const
    {
      return *Azimuthal_Mode_Number_pt;
    }

    /// Pointer to azimuthal mode number k in e^ik(theta) decomposition
    int*& azimuthal_mode_number_pt()
    {
      return Azimuthal_Mode_Number_pt;
    }

    /// Viscosity ratio for element: element's viscosity relative
    /// to the viscosity used in the definition of the Reynolds number
    const double& viscosity_ratio() const
    {
      return *Viscosity_Ratio_pt;
    }

    /// Pointer to the viscosity ratio
    double*& viscosity_ratio_pt()
    {
      return Viscosity_Ratio_pt;
    }

    /// Density ratio for element: element's density relative
    /// to the viscosity used in the definition of the Reynolds number
    const double& density_ratio() const
    {
      return *Density_Ratio_pt;
    }

    /// Pointer to the density ratio
    double*& density_ratio_pt()
    {
      return Density_Ratio_pt;
    }

    /// Access function for the base flow velocity pointer
    void (*&base_flow_u_fct_pt())(const double& time,
                                  const Vector<double>& x,
                                  Vector<double>& f)
    {
      return Base_flow_u_fct_pt;
    }

    /// Access function for the derivatives of the base flow
    /// w.r.t. global coordinates solution pointer
    void (*&base_flow_dudx_fct_pt())(const double& time,
                                     const Vector<double>& x,
                                     DenseMatrix<double>& f)
    {
      return Base_flow_dudx_fct_pt;
    }

    /// Access function for the derivatives of the base flow
    /// velocities w.r.t. time solution pointer
    void (*&base_flow_dudt_fct_pt())(const double& time,
                                     const Vector<double>& x,
                                     Vector<double>& f)
    {
      return Base_flow_dudt_fct_pt;
    }

    /// Access function for the base flow pressure pointer
    void (*&base_flow_p_fct_pt())(const double& time,
                                  const Vector<double>& x,
                                  double& f)
    {
      return Base_flow_p_fct_pt;
    }

    /// Access function for the derivs w.r.t. nodal coords X_{pq}
    /// of the spatial derivatives of base flow velocities pointer
    void (*&base_flow_d_dudx_dX_fct_pt())(const double& time,
                                          const Vector<double>& x,
                                          RankFourTensor<double>& f)
    {
      return Base_flow_d_dudx_dX_fct_pt;
    }

    /// Access function for the body-force pointer
    void (*&body_force_fct_pt())(const double& time,
                                 const Vector<double>& x,
                                 Vector<double>& f)
    {
      return Body_force_fct_pt;
    }

    /// Access function for the source-function pointer
    double (*&source_fct_pt())(const double& time, const Vector<double>& x)
    {
      return Source_fct_pt;
    }

    /// Access function for the derivatives of the base flow
    /// velocities w.r.t. local coordinates solution pointer
    void (*&base_flow_duds_fct_pt())(const double& time,
                                     const Vector<double>& x,
                                     DenseMatrix<double>& f)
    {
      return Base_flow_duds_fct_pt;
    }

    /// Return the number of pressure degrees of freedom
    /// associated with a single pressure component in the element
    virtual unsigned npres_lin_axi_nst() const = 0;

    /// Return the index at which the i-th component of the unknown
    /// perturbation to the nodal position is stored. The default value, i,
    /// is appropriate for single-physics problems. In derived multi-physics
    /// elements, this function should be overloaded to reflect the chosen
    /// storage scheme. Note that these equations require that the unknowns
    /// are always stored at the same indices at each node.
    virtual inline unsigned xhat_index_lin_axi_nst(const unsigned& n,
                                                   const unsigned& i) const
    {
      return i + 6;
    }

    /// Return the index at which the i-th unknown velocity
    /// component is stored. The default value, i+4, is appropriate for
    /// single-physics problems, since it stores the velocity components
    /// immediately after the four perturbations to the nodal positions.
    /// In derived multi-physics elements, this function should be
    /// overloaded to reflect the chosen storage scheme.
    /// Note that these equations require that the unknowns are always
    /// stored at the same indices at each node.
    virtual inline unsigned u_index_lin_axi_nst(const unsigned& i) const
    {
      return i;
    }

    /// Return the i-th component of du/dt at local node n.
    /// Uses suitably interpolated value for hanging nodes.
    double du_dt_lin_axi_nst(const unsigned& n, const unsigned& i) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt = this->node_pt(n)->time_stepper_pt();

      // Initialise dudt
      double dudt = 0.0;

      // Loop over the timesteps, if there is a non-steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Get the index at which the velocity is stored
        const unsigned u_nodal_index = u_index_lin_axi_nst(i);

        // Determine number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();

        // Add the contributions to the time derivative
        for (unsigned t = 0; t < n_time; t++)
        {
          dudt +=
            time_stepper_pt->weight(1, t) * nodal_value(t, n, u_nodal_index);
        }
      }

      return dudt;
    }

    /// Return the i-th component of dnodal_xhat/dt at local node n.
    /// Uses suitably interpolated value for hanging nodes.
    double dnodal_position_perturbation_dt_lin_axi_nst(const unsigned& n,
                                                       const unsigned& i) const
    {
      // Get the node's positional timestepper
      TimeStepper* time_stepper_pt =
        this->node_pt(n)->position_time_stepper_pt();

      // Initialise dxhat/dt
      double dnodal_xhat_dt = 0.0;

      // Loop over the timesteps, if there is a non-steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Get the index at which the perturbation to the nodal position
        // is stored
        const unsigned xhat_nodal_index = xhat_index_lin_axi_nst(n, i);

        // Determine number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();

        // Add the contributions to the time derivative
        for (unsigned t = 0; t < n_time; t++)
        {
          dnodal_xhat_dt +=
            time_stepper_pt->weight(1, t) * nodal_value(t, n, xhat_nodal_index);
        }
      }

      return dnodal_xhat_dt;
    }

    /// Disable ALE, i.e. assert the mesh is not moving -- you do this
    /// at your own risk!
    void disable_ALE()
    {
      ALE_is_disabled = true;
    }

    /// (Re-)enable ALE, i.e. take possible mesh motion into account
    /// when evaluating the time-derivative. Note: By default, ALE is
    /// enabled, at the expense of possibly creating unnecessary work
    /// in problems where the mesh is, in fact, stationary.
    void enable_ALE()
    {
      ALE_is_disabled = false;
    }

    /// Return the i-th pressure value at local pressure "node" n_p.
    /// Uses suitably interpolated value for hanging nodes.
    virtual double p_lin_axi_nst(const unsigned& n_p,
                                 const unsigned& i) const = 0;

    // Pin the pressure. i = 0 - cosine, i = 1 - sine
    // virtual void pin_pressure(const unsigned& i);

    /// Which nodal value represents the pressure?
    //  N.B. This function has return type "int" (rather than "unsigned"
    //  as in the u_index case) so that we can return the "magic" number
    //  "Pressure_not_stored_at_node" ( = -100 )
    virtual inline int p_index_lin_axi_nst(const unsigned& i) const
    {
      return Pressure_not_stored_at_node;
    }

    /// Get integral of kinetic energy over element plus deriv w.r.t. time
    void dkin_energy_dt(double& dkin_en_dt, double& kin_en) const;

    /// Strain-rate tensor: \f$ e_{ij} \f$
    /// where \f$ i,j = r,z,\theta \f$ (in that order)
    void strain_rate(const Vector<double>& s,
                     DenseMatrix<double>& strain_rate) const;

    /// Output function in tecplot format:
    /// r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S, R^C, R^S, Z^C, Z^S.
    /// Default number of plot points
    void output(std::ostream& outfile)
    {
      const unsigned nplot = 5;
      output(outfile, nplot);
    }

    /// Output function in tecplot format:
    /// r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S, R^C, R^S, Z^C, Z^S.
    /// Use nplot points in each coordinate direction
    void output(std::ostream& outfile, const unsigned& nplot);

    /// Output function in tecplot format:
    /// r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S, R^C, R^S, Z^C, Z^S.
    /// Default number of plot points
    void output(FILE* file_pt)
    {
      const unsigned nplot = 5;
      output(file_pt, nplot);
    }

    /// Output function in tecplot format:
    /// r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S, R^C, R^S, Z^C, Z^S.
    /// Use nplot points in each coordinate direction
    void output(FILE* file_pt, const unsigned& nplot);

    /// Output function: r, z, U^C, U^S, V^C, V^S, W^C, W^S,
    /// in tecplot format. nplot points in each coordinate direction
    /// at timestep t (t=0: present; t>0: previous timestep)
    void output_veloc(std::ostream& outfile,
                      const unsigned& nplot,
                      const unsigned& t);

    /// Compute the element's residual Vector
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // and using a dummy matrix argument
      fill_in_generic_residual_contribution_lin_axi_nst(
        residuals,
        GeneralisedElement::Dummy_matrix,
        GeneralisedElement::Dummy_matrix,
        0);
    }

    /// Compute the element's residual Vector and the jacobian matrix.
    /// Virtual function can be overloaded by hanging-node version.
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      // fill_in_generic_residual_contribution_lin_axi_nst(
      //   residuals, jacobian, GeneralisedElement::Dummy_matrix, 1);
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      this->fill_in_contribution_to_jacobian(residuals, jacobian);
      const unsigned n_dofs = this->ndof();
      Vector<double> dummy_residuals(n_dofs);
      DenseMatrix<double> dummy_jacobian(n_dofs, n_dofs);
      // Call the generic routine with the flag set to 2
      fill_in_generic_residual_contribution_lin_axi_nst(
        dummy_residuals, dummy_jacobian, mass_matrix, 2);
    }

    /// Return the i-th component of the FE interpolated velocity
    /// u[i] at local coordinate s
    double interpolated_u_lin_axi_nst(const Vector<double>& s,
                                      const unsigned& i) const
    {
      // Determine number of nodes in the element
      const unsigned n_node = nnode();

      // Provide storage for local shape functions
      Shape psi(n_node);

      // Find values of shape functions
      shape(s, psi);

      // Get the index at which the velocity is stored
      const unsigned u_nodal_index = u_index_lin_axi_nst(i);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += nodal_value(l, u_nodal_index) * psi[l];
      }

      return (interpolated_u);
    }

    /// Return the i-th component of the FE interpolated pressure
    /// p[i] at local coordinate s
    double interpolated_p_lin_axi_nst(const Vector<double>& s,
                                      const unsigned& i) const
    {
      // Determine number of pressure nodes in the element
      const unsigned n_pressure_nodes = npres_lin_axi_nst();

      // Provide storage for local shape functions
      Shape psi(n_pressure_nodes);

      // Find values of shape functions
      pshape_lin_axi_nst(s, psi);

      // Initialise value of p
      double interpolated_p = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_pressure_nodes; l++)
      {
        // N.B. The pure virtual function p_lin_axi_nst(...)
        // automatically calculates the index at which the pressure value
        // is stored, so we don't need to worry about this here
        interpolated_p += p_lin_axi_nst(l, i) * psi[l];
      }

      return (interpolated_p);
    }

    /// Return the i-th component of the FE interpolated
    /// perturbation to the nodal position xhat[i] at local coordinate s
    double interpolated_nodal_position_perturbation_lin_axi_nst(
      const Vector<double>& s, const unsigned& i) const
    {
      // Determine number of nodes in the element
      const unsigned n_node = nnode();

      // Provide storage for local shape functions
      Shape psi(n_node);

      // Find values of shape functions
      shape(s, psi);


      // Initialise value of xhat
      double interpolated_xhat = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the index at which the perturbation to the nodal position
        // is stored
        const unsigned xhat_nodal_index = xhat_index_lin_axi_nst(l, i);
        interpolated_xhat += nodal_value(l, xhat_nodal_index) * psi[l];
      }

      return (interpolated_xhat);
    }
  }; // End of FullLinearisedAxisymmetricNavierStokesEquations class definition

} // namespace oomph

#endif