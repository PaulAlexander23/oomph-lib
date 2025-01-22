#ifndef FULL_CONTINUATION_PROBLEM_HEADER
#define FULL_CONTINUATION_PROBLEM_HEADER

#include "singular_axisym_dynamic_cap_problem.h"
#include "height_continuation_problem.h"
#include "parameters.h"

namespace oomph
{
  template<class ELEMENT, class TIMESTEPPER>
  class FullContinuationProblem
    : public virtual SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>,
      public virtual HeightContinuationProblem
  {
  public:
    FullContinuationProblem(Params& params)
      : SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>(params),
        HeightContinuationProblem()
    {
    }

    void set_continuation_parameter(double& parameter)
    {
      this->set_continuation_parameter(parameter);
    }
  };
}; // namespace oomph

#endif
