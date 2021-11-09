// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//           Version 0.90. August 3, 2009.
// LIC//
// LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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

#ifndef PROJECTABLE_HELE_SHAW_ELEMENT_HEADER
#define PROJECTABLE_HELE_SHAW_ELEMENT_HEADER

// OOMPH-LIB headers
#include "../generic/Qelements.h"
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/projection.h"

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

namespace oomph
{
  //==========================================================
  /// HeleShaw upgraded to become projectable
  //==========================================================
  template<class HELE_SHAW_ELEMENT>
  class ProjectableHeleShawElement
    : public virtual ProjectableElement<HELE_SHAW_ELEMENT>
  {
  public:
    /// \short Specify the values associated with field fld.
    /// The information is returned in a vector of pairs which comprise
    /// the Data object and the value within it, that correspond to field fld.
    Vector<std::pair<Data*, unsigned>> data_values_of_field(
      const unsigned& fld);

    /// \short Number of fields to be projected: Just one
    unsigned nfields_for_projection();

    /// \short Number of history values to be stored for fld-th field.
    unsigned nhistory_values_for_projection(const unsigned& fld);

    ///\short Number of positional history values
    unsigned nhistory_values_for_coordinate_projection();

    /// \short Return Jacobian of mapping and shape functions of field fld
    /// at local coordinate s
    double jacobian_and_shape_of_field(const unsigned& fld,
                                       const Vector<double>& s,
                                       Shape& psi);

    /// \short Return interpolated field fld at local coordinate s, at time
    /// level t (t=0: present; t>0: history values)
    double get_field(const unsigned& t,
                     const unsigned& fld,
                     const Vector<double>& s);

    /// Return number of values in field fld: One per node
    unsigned nvalue_of_field(const unsigned& fld);

    /// Return local equation number of value j in field fld.
    int local_equation(const unsigned& fld, const unsigned& j);
  };

} // namespace oomph

#include "projectable_hele_shaw_element.template.cc"

#endif
