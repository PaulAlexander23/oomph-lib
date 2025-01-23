#ifndef FULL_CONTINUATION_PROBLEM_HEADER
#define FULL_CONTINUATION_PROBLEM_HEADER

#include "singular_axisym_dynamic_cap_problem.h"
#include "parameters.h"

namespace oomph
{
  template<class ELEMENT, class TIMESTEPPER>
  class FullContinuationProblem
    : public virtual SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>
  {
  public:
    FullContinuationProblem(Params* const& params)
      : SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>(params)
    {
    }
  };
}; // namespace oomph

#endif
