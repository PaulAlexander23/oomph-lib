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
  private:
    /// Storage for the height element mesh
    Mesh* Height_mesh_pt;

    /// Storage for the traded height step data
    Data* Traded_data_pt;

  public:
    FullContinuationProblem(Params* const& params)
      : SingularAxisymDynamicCapProblem<ELEMENT, TIMESTEPPER>(params),
        Height_mesh_pt(new Mesh),
        Traded_data_pt(new Data(1))
    {
      this->add_sub_mesh(Height_mesh_pt);
      this->add_global_data(Traded_data_pt);
      this->rebuild_global_mesh();

      Traded_data_pt->pin(0);
    }
  };
}; // namespace oomph

#endif
