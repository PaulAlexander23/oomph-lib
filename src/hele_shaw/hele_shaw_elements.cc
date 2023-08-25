#include "hele_shaw_elements.h"

#ifndef OOMPH_HELE_SHAW_ELEMENTS_TEMPLATE
#define OOMPH_HELE_SHAW_ELEMENTS_TEMPLATE

namespace oomph
{
  //======================================================================
  /// QHeleShawElement elements are linear/quadrilateral/brick-shaped
  /// HeleShaw elements with isoparametric interpolation for the function.
  //======================================================================
  // JACK - SPECIFIC ELEMENT

  /// \short Output function:
  ///  x,y,u   or    x,y,z,u
  template<unsigned NNODE_1D>
  void QHeleShawElement<NNODE_1D>::output(std::ostream& outfile)
  {
    HeleShawEquations::output(outfile);
  }

  ///  \short Output function:
  ///   x,y,u   or    x,y,z,u at n_plot^2 plot points
  template<unsigned NNODE_1D>
  void QHeleShawElement<NNODE_1D>::output(std::ostream& outfile,
                                          const unsigned& n_plot)
  {
    HeleShawEquations::output(outfile, n_plot);
  }

  /// \short C-style output function:
  ///  x,y,u   or    x,y,z,u
  template<unsigned NNODE_1D>
  void QHeleShawElement<NNODE_1D>::output(FILE* file_pt)
  {
    HeleShawEquations::output(file_pt);
  }

  ///  \short C-style output function:
  ///   x,y,u   or    x,y,z,u at n_plot^2 plot points
  template<unsigned NNODE_1D>
  void QHeleShawElement<NNODE_1D>::output(FILE* file_pt, const unsigned& n_plot)
  {
    HeleShawEquations::output(file_pt, n_plot);
  }

  /// \short Output function for an exact solution:
  ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^2 plot points
  template<unsigned NNODE_1D>
  void QHeleShawElement<NNODE_1D>::output_fct(
    std::ostream& outfile,
    const unsigned& n_plot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    HeleShawEquations::output_fct(outfile, n_plot, exact_soln_pt);
  }

  /// \short Output function for a time-dependent exact solution.
  ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^2 plot points
  /// (Calls the steady version)
  template<unsigned NNODE_1D>
  void QHeleShawElement<NNODE_1D>::output_fct(
    std::ostream& outfile,
    const unsigned& n_plot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
    HeleShawEquations::output_fct(outfile, n_plot, time, exact_soln_pt);
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //======================================================================
  /// Set the data for the number of Variables at each node, always one
  /// in every case
  //======================================================================
  template<unsigned NNODE_1D>
  const unsigned QHeleShawElement<NNODE_1D>::Initial_Nvalue = 1;


  //====================================================================
  /// Force build of templates
  //====================================================================
  template class QHeleShawElement<2>;
  template class QHeleShawElement<3>;
  template class QHeleShawElement<4>;

} // namespace oomph
#endif
