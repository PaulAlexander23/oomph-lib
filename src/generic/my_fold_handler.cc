#include "my_fold_handler.h"
#include "elements.h"
#include "problem.h"
#include "mesh.h"

namespace oomph
{
  //========================================================================
  /// Constructor: Initialise the fold handler,
  /// by setting initial guesses for Y, Phi and calculating Count.
  /// If the system changes, a new  handler must be constructed.
  //========================================================================
  MyFoldHandler::MyFoldHandler(Problem* const& problem_pt,
                               double* const& parameter_pt,
                               double* const& eigenvalue_pt)
    : Solve_which_system(Full_augmented),
      Parameter_pt(parameter_pt),
      Eigenvalue_pt(eigenvalue_pt)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of degrees of freedom
    Ndof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution* dist_pt =
      new LinearAlgebraDistribution(problem_pt->communicator_pt(), Ndof, false);

    // Resize the vectors of additional dofs and constants
    Phi.resize(Ndof);
    Y.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in the element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Calculate the value Phi by
    // solving the system JPhi = dF/dlambda

    // Locally cache the linear solver
    LinearSolver* const linear_solver_pt = problem_pt->linear_solver_pt();

    // Save the status before entry to this routine
    bool enable_resolve = linear_solver_pt->is_resolve_enabled();

    // We need to do a resolve
    linear_solver_pt->enable_resolve();

    // Storage for the solution
    DoubleVector x(dist_pt, 0.0);

    // Solve the standard problem, we only want to make sure that
    // we factorise the matrix, if it has not been factorised. We shall
    // ignore the return value of x.
    linear_solver_pt->solve(problem_pt, x);

    // Get the vector dresiduals/dparameter
    problem_pt->get_derivative_wrt_global_parameter(parameter_pt, x);

    // Copy rhs vector into local storage so it doesn't get overwritten
    // if the linear solver decides to initialise the solution vector, say,
    // which it's quite entitled to do!
    DoubleVector input_x(x);

    // Now resolve the system with the new RHS and overwrite the solution
    linear_solver_pt->resolve(input_x, x);

    // Restore the storage status of the linear solver
    if (enable_resolve)
    {
      linear_solver_pt->enable_resolve();
    }
    else
    {
      linear_solver_pt->disable_resolve();
    }

    // Add the global parameter as an unknown to the problem
    problem_pt->Dof_pt.push_back(parameter_pt);


    // Normalise the initial guesses for phi
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += x[n] * x[n];
    }
    length = sqrt(length);

    // Now add the null space components to the problem unknowns
    // and initialise them and Phi to the same normalised values
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Y[n]);
      Y[n] = Phi[n] = -x[n] / length;
    }

    // delete the dist_pt
    problem_pt->Dof_distribution_pt->build(
      problem_pt->communicator_pt(), Ndof * 2 + 1, true);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

    delete dist_pt;
  }


  /// Constructor in which the eigenvector is passed as an initial
  /// guess
  MyFoldHandler::MyFoldHandler(Problem* const& problem_pt,
                               double* const& parameter_pt,
                               double* const& eigenvalue_pt,
                               const DoubleVector& eigenvector)
    : Solve_which_system(Full_augmented),
      Parameter_pt(parameter_pt),
      Eigenvalue_pt(eigenvalue_pt)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of degrees of freedom
    Ndof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution* dist_pt =
      new LinearAlgebraDistribution(problem_pt->communicator_pt(), Ndof, false);

    // Resize the vectors of additional dofs and constants
    Phi.resize(Ndof);
    Y.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in the element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Add the global parameter as an unknown to the problem
    problem_pt->Dof_pt.push_back(parameter_pt);


    // Normalise the initial guesses for the eigenvecor
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += eigenvector[n] * eigenvector[n];
    }
    length = sqrt(length);

    // Now add the null space components to the problem unknowns
    // and initialise them and Phi to the same normalised values
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Y[n]);
      Y[n] = Phi[n] = eigenvector[n] / length;
    }

    // delete the dist_pt
    problem_pt->Dof_distribution_pt->build(
      problem_pt->communicator_pt(), Ndof * 2 + 1, true);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

    delete dist_pt;
  }


  /// Constructor in which the eigenvector and normalisation
  /// vector are  passed as an initial guess
  MyFoldHandler::MyFoldHandler(Problem* const& problem_pt,
                               double* const& parameter_pt,
                               double* const& eigenvalue_pt,
                               const DoubleVector& eigenvector,
                               const DoubleVector& normalisation)
    : Solve_which_system(Full_augmented),
      Parameter_pt(parameter_pt),
      Eigenvalue_pt(eigenvalue_pt)
  {
    // Set the problem pointer
    Problem_pt = problem_pt;
    // Set the number of degrees of freedom
    Ndof = problem_pt->ndof();

    // create the linear algebra distribution for this solver
    // currently only global (non-distributed) distributions are allowed
    LinearAlgebraDistribution* dist_pt =
      new LinearAlgebraDistribution(problem_pt->communicator_pt(), Ndof, false);

    // Resize the vectors of additional dofs and constants
    Phi.resize(Ndof);
    Y.resize(Ndof);
    Count.resize(Ndof, 0);

    // Loop over all the elements in the problem
    unsigned n_element = problem_pt->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      GeneralisedElement* elem_pt = problem_pt->mesh_pt()->element_pt(e);
      // Loop over the local freedoms in the element
      unsigned n_var = elem_pt->ndof();
      for (unsigned n = 0; n < n_var; n++)
      {
        // Increase the associated global equation number counter
        ++Count[elem_pt->eqn_number(n)];
      }
    }

    // Add the global parameter as an unknown to the problem
    problem_pt->Dof_pt.push_back(parameter_pt);


    // Normalise the initial guesses for the eigenvecor
    double length = 0.0;
    for (unsigned n = 0; n < Ndof; n++)
    {
      length += eigenvector[n] * normalisation[n];
    }
    length = sqrt(length);

    // Now add the null space components to the problem unknowns
    // and initialise them and Phi to the same normalised values
    for (unsigned n = 0; n < Ndof; n++)
    {
      problem_pt->Dof_pt.push_back(&Y[n]);
      Y[n] = eigenvector[n] / length;
      Phi[n] = normalisation[n];
    }

    // delete the dist_pt
    problem_pt->Dof_distribution_pt->build(
      problem_pt->communicator_pt(), Ndof * 2 + 1, true);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

    delete dist_pt;
  }


  //=================================================================
  /// The number of unknowns is 2n+1
  //================================================================
  unsigned MyFoldHandler::ndof(GeneralisedElement* const& elem_pt)
  {
    unsigned raw_ndof = elem_pt->ndof();
    // Return different values depending on the type of block decomposition
    switch (Solve_which_system)
    {
      case Full_augmented:
        return (2 * raw_ndof + 1);
        break;

      case Block_augmented_J:
        return (raw_ndof + 1);
        break;

      case Block_J:
        return raw_ndof;
        break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }

  //=====================================================================
  /// Return the global equation number associated with local equation
  /// number ieqn_local. We choose to number the unknowns according
  /// to the augmented system.
  //=======================================================================
  unsigned long MyFoldHandler::eqn_number(GeneralisedElement* const& elem_pt,
                                          const unsigned& ieqn_local)
  {
    // Find the number of non-augmented dofs in the element
    unsigned raw_ndof = elem_pt->ndof();
    // Storage for the global eqn number
    unsigned long global_eqn = 0;
    // If we are a "standard" unknown, just return the standard equation number
    if (ieqn_local < raw_ndof)
    {
      global_eqn = elem_pt->eqn_number(ieqn_local);
    }
    // Otherwise if we are at an unknown corresponding to the bifurcation
    // parameter return the global equation number of the parameter
    else if (ieqn_local == raw_ndof)
    {
      global_eqn = Ndof;
    }
    // Otherwise we are in the unknown corresponding to a null vector
    // return the global unknown Ndof + 1 + global unknown of "standard"
    // unknown.
    else
    {
      global_eqn = Ndof + 1 + elem_pt->eqn_number(ieqn_local - 1 - raw_ndof);
    }

    // Return the global equation number
    return global_eqn;
  }

  //====================================================================
  /// Formulate the augmented system
  //====================================================================
  void MyFoldHandler::get_residuals(GeneralisedElement* const& elem_pt,
                                    Vector<double>& residuals)
  {
    // Need to get raw residuals and jacobian
    unsigned raw_ndof = elem_pt->ndof();

    // Find out which system we are solving
    switch (Solve_which_system)
    {
        // If we are solving the standard system
      case Block_J:
      {
        // Get the basic residuals
        elem_pt->get_residuals(residuals);
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the basic residuals
        elem_pt->get_residuals(residuals);

        // Zero the final residual
        residuals[raw_ndof] = 0.0;
      }
      break;

      // If we are solving the full augmented system
      case Full_augmented:
      {
        DenseMatrix<double> jacobian(raw_ndof);
        // Get the basic residuals and jacobian initially
        elem_pt->get_jacobian(residuals, jacobian);

        // The normalisation equation must be initialised to -1.0/number of
        // elements
        residuals[raw_ndof] = -1.0 / Problem_pt->mesh_pt()->nelement();

        // Now assemble the equations Jy = cy
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          residuals[raw_ndof + 1 + i] =
            -(*Eigenvalue_pt) * Y[elem_pt->eqn_number(i)];
          // std::cout << *Eigenvalue_pt << std::endl;
          for (unsigned j = 0; j < raw_ndof; j++)
          {
            residuals[raw_ndof + 1 + i] +=
              jacobian(i, j) * Y[elem_pt->eqn_number(j)];
          }
          // Add the contribution to phi.y=1
          // Need to divide by the number of elements that contribute to this
          // unknown so that we do get phi.y exactly.
          unsigned global_eqn = elem_pt->eqn_number(i);
          residuals[raw_ndof] +=
            (Phi[global_eqn] * Y[global_eqn]) / Count[global_eqn];
        }
        std::cout << residuals[raw_ndof] << std::endl;
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=================================================================
  /// Calculate the elemental Jacobian matrix "d equation
  /// / d variable" for the augmented system
  //=================================================================
  void MyFoldHandler::get_jacobian(GeneralisedElement* const& elem_pt,
                                   Vector<double>& residuals,
                                   DenseMatrix<double>& jacobian)
  {
    // Find the number of augmented dofs
    unsigned augmented_ndof = ndof(elem_pt);
    // Find the non-augmented dofs
    unsigned raw_ndof = elem_pt->ndof();

    // Which system are we solving
    switch (Solve_which_system)
    {
        // If we are solving the original system
      case Block_J:
      {
        // Just get the raw jacobian and residuals
        elem_pt->get_jacobian(residuals, jacobian);
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the full residuals, we need them
        get_residuals(elem_pt, residuals);

        // Need to get the raw jacobian (and raw residuals)
        Vector<double> newres(augmented_ndof);
        elem_pt->get_jacobian(newres, jacobian);

        // Now do finite differencing stuff
        const double FD_step = 1.0e-8;
        // Fill in the first lot of finite differences
        {
          // increase the global parameter
          double* unknown_pt = Problem_pt->Dof_pt[Ndof];
          double init = *unknown_pt;
          *unknown_pt += FD_step;

          // Now do any possible updates
          Problem_pt->actions_after_change_in_bifurcation_parameter();

          // Get the new (modified) residuals
          get_residuals(elem_pt, newres);

          // The final column  is given by the difference
          // between the residuals
          for (unsigned n = 0; n < raw_ndof; n++)
          {
            jacobian(n, augmented_ndof - 1) =
              (newres[n] - residuals[n]) / FD_step;
          }
          // Reset the global parameter
          *unknown_pt = init;

          // Now do any possible updates
          Problem_pt->actions_after_change_in_bifurcation_parameter();
        }

        // Fill in the bottom row
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          unsigned local_eqn = elem_pt->eqn_number(n);
          jacobian(augmented_ndof - 1, n) = Phi[local_eqn] / Count[local_eqn];
        }
      }
      break;

        // Otherwise solving the full system
      case Full_augmented:
      {
        // Get the residuals for the augmented system
        get_residuals(elem_pt, residuals);

        // Need to get the raw residuals and jacobian
        Vector<double> newres(raw_ndof);
        DenseMatrix<double> newjac(raw_ndof);
        elem_pt->get_jacobian(newres, jacobian);

        // Fill in the jacobian on the diagonal sub-block of
        // the null-space equations
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          for (unsigned m = 0; m < raw_ndof; m++)
          {
            jacobian(raw_ndof + 1 + n, raw_ndof + 1 + m) = jacobian(n, m);
          }
        }

        // Now finite difference wrt the global unknown
        const double FD_step = 1.0e-8;
        // Fill in the first lot of finite differences
        {
          // increase the global parameter
          double* unknown_pt = Problem_pt->Dof_pt[Ndof];
          double init = *unknown_pt;
          *unknown_pt += FD_step;
          // Need to update the function
          Problem_pt->actions_after_change_in_bifurcation_parameter();

          // Get the new raw residuals and jacobian
          elem_pt->get_jacobian(newres, newjac);

          // The end of the first row is given by the difference
          // between the residuals
          for (unsigned n = 0; n < raw_ndof; n++)
          {
            jacobian(n, raw_ndof) = (newres[n] - residuals[n]) / FD_step;
            // The end of the second row is given by the difference multiplied
            // by the product
            for (unsigned l = 0; l < raw_ndof; l++)
            {
              jacobian(raw_ndof + 1 + n, raw_ndof) +=
                (newjac(n, l) - jacobian(n, l)) * Y[elem_pt->eqn_number(l)] /
                FD_step;
            }
          }
          // Reset the global parameter
          *unknown_pt = init;

          // Need to update the function
          Problem_pt->actions_after_change_in_bifurcation_parameter();
        }

        // Now fill in the first column of the second rows
        {
          for (unsigned n = 0; n < raw_ndof; n++)
          {
            unsigned long global_eqn = eqn_number(elem_pt, n);
            // Increase the first lot
            double* unknown_pt = Problem_pt->Dof_pt[global_eqn];
            double init = *unknown_pt;
            *unknown_pt += FD_step;
            Problem_pt->actions_before_newton_convergence_check(); /// ALICE

            // Get the new jacobian
            elem_pt->get_jacobian(newres, newjac);

            // Work out the differences
            for (unsigned k = 0; k < raw_ndof; k++)
            {
              // jacobian(raw_ndof+k,n) = 0.0;
              for (unsigned l = 0; l < raw_ndof; l++)
              {
                jacobian(raw_ndof + 1 + k, n) +=
                  (newjac(k, l) - jacobian(k, l)) * Y[elem_pt->eqn_number(l)] /
                  FD_step;
              }
            }
            *unknown_pt = init;
            Problem_pt->actions_before_newton_convergence_check(); /// ALICE
          }
        }

        // Fill in the row corresponding to the parameter
        for (unsigned n = 0; n < raw_ndof; n++)
        {
          unsigned global_eqn = elem_pt->eqn_number(n);
          jacobian(raw_ndof, raw_ndof + 1 + n) =
            Phi[global_eqn] / Count[global_eqn];
        }

        //   //Loop over the dofs
        //     for(unsigned n=0;n<augmented_ndof;n++)
        //      {
        //       unsigned long global_eqn = eqn_number(elem_pt,n);
        //       double* unknown_pt = Problem_pt->Dof_pt[global_eqn];
        //       double init = *unknown_pt;
        //       *unknown_pt += FD_step;

        //       //Get the new residuals
        //       get_residuals(elem_pt,newres);

        //       for(unsigned m=0;m<augmented_ndof;m++)
        //        {
        //         jacobian(m,n) = (newres[m] - residuals[m])/FD_step;
        //        }
        //       //Reset the unknown
        //       *unknown_pt = init;
        //      }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //====================================================================
  /// Formulate the derivatives of the augmented system with respect
  /// to a parameter
  //====================================================================
  void MyFoldHandler::get_dresiduals_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam)
  {
    // Need to get raw residuals and jacobian
    unsigned raw_ndof = elem_pt->ndof();

    // Find out which system we are solving
    switch (Solve_which_system)
    {
        // If we are solving the standard system
      case Block_J:
      {
        // Get the basic residual derivatives
        elem_pt->get_dresiduals_dparameter(parameter_pt, dres_dparam);
      }
      break;

      // If we are solving the augmented-by-one system
      case Block_augmented_J:
      {
        // Get the basic residual derivatives
        elem_pt->get_dresiduals_dparameter(parameter_pt, dres_dparam);

        // Zero the final derivative
        dres_dparam[raw_ndof] = 0.0;
      }
      break;

      // If we are solving the full augmented system
      case Full_augmented:
      {
        DenseMatrix<double> djac_dparam(raw_ndof);
        // Get the basic residuals and jacobian derivatives initially
        elem_pt->get_djacobian_dparameter(
          parameter_pt, dres_dparam, djac_dparam);

        // The normalisation equation does not depend on the parameter
        dres_dparam[raw_ndof] = 0.0;

        // Now assemble the equations dJy/dparameter = 0
        for (unsigned i = 0; i < raw_ndof; i++)
        {
          dres_dparam[raw_ndof + 1 + i] = 0.0;
          for (unsigned j = 0; j < raw_ndof; j++)
          {
            dres_dparam[raw_ndof + 1 + i] +=
              djac_dparam(i, j) * Y[elem_pt->eqn_number(j)];
          }
        }
      }
      break;

      default:
        std::ostringstream error_stream;
        error_stream
          << "The Solve_which_system flag can only take values 0, 1, 2"
          << " not " << Solve_which_system << "\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //========================================================================
  /// Overload the derivative of the residuals and jacobian
  /// with respect to a parameter so that it breaks because it should not
  /// be required
  //========================================================================
  void MyFoldHandler::get_djacobian_dparameter(
    GeneralisedElement* const& elem_pt,
    double* const& parameter_pt,
    Vector<double>& dres_dparam,
    DenseMatrix<double>& djac_dparam)
  {
    std::ostringstream error_stream;
    error_stream
      << "This function has not been implemented because it is not required\n";
    error_stream << "in standard problems.\n";
    error_stream
      << "If you find that you need it, you will have to implement it!\n\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //=====================================================================
  /// Overload the hessian vector product function so that
  /// it breaks because it should not be required
  //========================================================================
  void MyFoldHandler::get_hessian_vector_products(
    GeneralisedElement* const& elem_pt,
    Vector<double> const& Y,
    DenseMatrix<double> const& C,
    DenseMatrix<double>& product)
  {
    std::ostringstream error_stream;
    error_stream
      << "This function has not been implemented because it is not required\n";
    error_stream << "in standard problems.\n";
    error_stream
      << "If you find that you need it, you will have to implement it!\n\n";

    throw OomphLibError(
      error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  //==========================================================================
  /// Return the eigenfunction(s) associated with the bifurcation that
  /// has been detected in bifurcation tracking problems
  //==========================================================================
  void MyFoldHandler::get_eigenfunction(Vector<DoubleVector>& eigenfunction)
  {
    // There is only one (real) null vector
    eigenfunction.resize(1);
    LinearAlgebraDistribution dist(Problem_pt->communicator_pt(), Ndof, false);
    // Rebuild the vector
    eigenfunction[0].build(&dist, 0.0);
    // Set the value to be the null vector
    for (unsigned n = 0; n < Ndof; n++)
    {
      eigenfunction[0][n] = Y[n];
    }
  }


  //=======================================================================
  /// Destructor return the problem to its original (non-augmented) state
  //=======================================================================
  MyFoldHandler::~MyFoldHandler()
  {
    // If we are using the block solver reset the problem's linear solver
    // to the original one
    AugmentedBlockFoldLinearSolver* block_fold_solver_pt =
      dynamic_cast<AugmentedBlockFoldLinearSolver*>(
        Problem_pt->linear_solver_pt());

    if (block_fold_solver_pt)
    {
      // Reset the problem's linear solver
      Problem_pt->linear_solver_pt() = block_fold_solver_pt->linear_solver_pt();
      // Delete the block solver
      delete block_fold_solver_pt;
    }

    // Resize the number of dofs
    Problem_pt->Dof_pt.resize(Ndof);
    Problem_pt->Dof_distribution_pt->build(
      Problem_pt->communicator_pt(), Ndof, false);
    // Remove all previous sparse storage used during Jacobian assembly
    Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
  }

  //====================================================================
  /// Set to solve the augmented-by-one block system.
  //===================================================================
  void MyFoldHandler::solve_augmented_block_system()
  {
    // Only bother to do anything if we haven't already set the flag
    if (Solve_which_system != Block_augmented_J)
    {
      // If we were solving the system with the original jacobian add the
      // parameter
      if (Solve_which_system == Block_J)
      {
        Problem_pt->Dof_pt.push_back(Parameter_pt);
      }

      // Restrict the problem to the standard variables and
      // the bifurcation parameter only
      Problem_pt->Dof_pt.resize(Ndof + 1);
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof + 1, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);
      // Set the solve flag
      Solve_which_system = Block_augmented_J;
    }
  }


  //====================================================================
  /// Set to solve the block system. The boolean flag specifies
  /// whether the block decomposition should use exactly the same jacobian
  //===================================================================
  void MyFoldHandler::solve_block_system()
  {
    // Only bother to do anything if we haven't already set the flag
    if (Solve_which_system != Block_J)
    {
      // Restrict the problem to the standard variables
      Problem_pt->Dof_pt.resize(Ndof);
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

      // Set the solve flag
      Solve_which_system = Block_J;
    }
  }

  //=================================================================
  /// Set to Solve non-block system
  //=================================================================
  void MyFoldHandler::solve_full_system()
  {
    // Only do something if we are not solving the full system
    if (Solve_which_system != Full_augmented)
    {
      // If we were solving the problem without any augmentation,
      // add the parameter again
      if (Solve_which_system == Block_J)
      {
        Problem_pt->Dof_pt.push_back(Parameter_pt);
      }

      // Always add the additional unknowns back into the problem
      for (unsigned n = 0; n < Ndof; n++)
      {
        Problem_pt->Dof_pt.push_back(&Y[n]);
      }

      // update the Dof distribution pt
      Problem_pt->Dof_distribution_pt->build(
        Problem_pt->communicator_pt(), Ndof * 2 + 1, false);
      // Remove all previous sparse storage used during Jacobian assembly
      Problem_pt->Sparse_assemble_with_arrays_previous_allocation.resize(0);

      Solve_which_system = Full_augmented;
    }
  }
} // namespace oomph
