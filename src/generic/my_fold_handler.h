#ifndef MY_FOLD_HANDLER_HEADER
#define MY_FOLD_HANDLER_HEADER

//#include "generic.h"
#include "matrices.h"
#include "linear_solver.h"
#include "double_vector_with_halo.h"
#include "assembly_handler.h"

namespace oomph
{
  // Forward class definition of the element
  class GeneralisedElement;

  // Forward class definition of the problem
  class Problem;

  //===================================================================
  /// A class that is used to assemble the augmented system that defines
  /// a fold (saddle-node) or limit point. The "standard" problem must
  /// be a function of a global paramter \f$\lambda\f$, and a
  /// solution is \f$R(u,\lambda) = 0 \f$ , where \f$ u \f$ are the unknowns
  /// in the problem. A limit point is formally specified by the augmented
  /// system of size \f$ 2N+1 \f$
  /// \f[ R(u,\lambda) = 0, \f]
  /// \f[ Jy = c y, \f]
  /// \f[\phi\cdot y = 1. \f]
  /// In the above \f$ J \f$ is the usual Jacobian matrix, \f$ dR/du \f$.
  /// \f$ c \f$ is a constant eigenvalue and
  /// \f$ \phi \f$ is a constant vector that is chosen to
  /// ensure that the null vector, \f$ y \f$, is not trivial.
  //======================================================================
  class MyFoldHandler : public AssemblyHandler
  {
    friend class AugmentedBlockFoldLinearSolver;

    /// A little private enum to determine whether we are solving
    /// the block system or not
    enum
    {
      Full_augmented,
      Block_J,
      Block_augmented_J
    };

    /// Integer flag to indicate which system should be assembled.
    /// There are three possibilities. The full augmented system (0),
    /// the non-augmented jacobian system (1), a system in which the
    /// jacobian is augmented by 1 row and column to ensure that it is
    /// non-singular (2). See the enum above
    unsigned Solve_which_system;

    /// Pointer to the problem
    Problem* Problem_pt;

    /// Store the number of degrees of freedom in the non-augmented
    /// problem
    unsigned Ndof;

    /// A constant vector used to ensure that the null vector
    /// is not trivial
    Vector<double> Phi;

    /// Storage for the null vector
    Vector<double> Y;

    /// A vector that is used to determine how many elements
    /// contribute to a particular equation. It is used to ensure
    /// that the global system is correctly formulated.
    Vector<int> Count;

    /// Storage for the pointer to the parameter
    double* Parameter_pt;

    double* Eigenvalue_pt;

  public:
    /// Constructor:
    /// initialise the fold handler, by setting initial guesses
    /// for Y, Phi and calculating count. If the system changes, a new
    /// fold handler must be constructed
    MyFoldHandler(Problem* const& problem_pt,
                  double* const& parameter_pt,
                  double* const& eigenvalue_pt);


    /// Constructor in which initial eigenvector can be passed
    MyFoldHandler(Problem* const& problem_pt,
                  double* const& parameter_pt,
                  double* const& eigenvalue_pt,
                  const DoubleVector& eigenvector);

    /// Constructor in which initial eigenvector
    /// and normalisation can be passed
    MyFoldHandler(Problem* const& problem_pt,
                  double* const& parameter_pt,
                  double* const& eigenvalue_pt,
                  const DoubleVector& eigenvector,
                  const DoubleVector& normalisation);


    /// Destructor, return the problem to its original state
    /// before the augmented system was added
    ~MyFoldHandler();

    /// Get the number of elemental degrees of freedom
    unsigned ndof(GeneralisedElement* const& elem_pt);

    /// Get the global equation number of the local unknown
    unsigned long eqn_number(GeneralisedElement* const& elem_pt,
                             const unsigned& ieqn_local);

    /// Get the residuals
    void get_residuals(GeneralisedElement* const& elem_pt,
                       Vector<double>& residuals);

    /// Calculate the elemental Jacobian matrix "d equation
    /// / d variable".
    void get_jacobian(GeneralisedElement* const& elem_pt,
                      Vector<double>& residuals,
                      DenseMatrix<double>& jacobian);

    /// Overload the derivatives of the residuals with respect to
    /// a parameter to apply to the augmented system
    void get_dresiduals_dparameter(GeneralisedElement* const& elem_pt,
                                   double* const& parameter_pt,
                                   Vector<double>& dres_dparam);

    /// Overload the derivative of the residuals and jacobian
    /// with respect to a parameter so that it breaks
    void get_djacobian_dparameter(GeneralisedElement* const& elem_pt,
                                  double* const& parameter_pt,
                                  Vector<double>& dres_dparam,
                                  DenseMatrix<double>& djac_dparam);

    /// Overload the hessian vector product function so that
    /// it breaks
    void get_hessian_vector_products(GeneralisedElement* const& elem_pt,
                                     Vector<double> const& Y,
                                     DenseMatrix<double> const& C,
                                     DenseMatrix<double>& product);

    /// Indicate that we are tracking a fold bifurcation by returning 1
    int bifurcation_type() const
    {
      return 1;
    }

    /// Return a pointer to the
    /// bifurcation parameter in bifurcation tracking problems
    double* bifurcation_parameter_pt() const
    {
      return Parameter_pt;
    }

    /// Return the eigenfunction(s) associated with the bifurcation that
    /// has been detected in bifurcation tracking problems
    void get_eigenfunction(Vector<DoubleVector>& eigenfunction);

    /// Set to solve the augmented block system
    void solve_augmented_block_system();

    /// Set to solve the block system
    void solve_block_system();

    /// Solve non-block system
    void solve_full_system();
  };
} // namespace oomph
#endif
