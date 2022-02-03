#include "generic.h"

namespace oomph
{
  // A Comparison operator for the boundary nodes
  class CompareNodeCoordinatesX
  {
  public:
    /// The actual comparison operator
    int operator()(Node* const& node1_pt, Node* const& node2_pt)
    {
      unsigned n_dim = node1_pt->ndim();
      if (n_dim != node2_pt->ndim())
      {
        throw OomphLibError("Can't compare two nodes of different dimension",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Make sure to handle the finite precision problems
      return node1_pt->x(0) < node2_pt->x(0);
    }
  };

  // A Comparison operator for the boundary nodes
  class CompareNodeCoordinatesY
  {
  public:
    /// The actual comparison operator
    int operator()(Node* const& node1_pt, Node* const& node2_pt)
    {
      unsigned n_dim = node1_pt->ndim();
      if (n_dim != node2_pt->ndim())
      {
        throw OomphLibError("Can't compare two nodes of different dimension",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Make sure to handle the finite precision problems
      return node1_pt->x(1) < node2_pt->x(1);
    }
  };
} // namespace oomph
