#ifndef ELEMENT_WITH_ERROR_HEADER
#define ELEMENT_WITH_ERROR_HEADER

#include "generic.h"

namespace oomph
{
  template<class ELEMENT>
  class ElementWithError : public virtual ELEMENT
  {
    double Error;

  public:
    ElementWithError() : ELEMENT(), Error(0.0) {}

    /// Set error value for post-processing
    void set_error(const double& error)
    {
      Error = error;
    }

    void get_error(double& error)
    {
      error = Error;
    }
  };

  template<class ELEMENT>
  class FaceGeometry<ElementWithError<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  };

  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<ElementWithError<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  };


}; // namespace oomph
#endif
