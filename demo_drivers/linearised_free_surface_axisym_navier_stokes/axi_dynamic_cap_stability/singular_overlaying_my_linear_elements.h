#ifndef SINGULAR_OVERLAYING_MY_LINEAR_ELEMENT_HEADER
#define SINGULAR_OVERLAYING_MY_LINEAR_ELEMENT_HEADER

#include "overlaying_my_linear_element.h"

namespace oomph
{
  template<class BASE_ELEMENT>
  class SingularOverlayingMyLinearElement
    : public virtual OverlayingMyLinearElement<BASE_ELEMENT>
  {
  private:
    bool IsAugmented;

  public:
    SingularOverlayingMyLinearElement()
      : OverlayingMyLinearElement<BASE_ELEMENT>(), IsAugmented(false)
    {
    }

    bool is_augmented() const
    {
      return IsAugmented;
    }

    void augment()
    {
      IsAugmented = true;
    }
  };

  template<class BASE_ELEMENT>
  class FaceGeometry<SingularOverlayingMyLinearElement<BASE_ELEMENT>>
    : public TElement<1, 3>
  {
  };
}; // namespace oomph

#endif
