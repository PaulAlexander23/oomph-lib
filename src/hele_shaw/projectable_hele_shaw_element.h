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

  //==========================================================
  /// HeleShaw upgraded to become projectable
  //==========================================================
  template<class HELE_SHAW_ELEMENT>
  Vector<std::pair<Data*, unsigned>> ProjectableHeleShawElement<
    HELE_SHAW_ELEMENT>::data_values_of_field(const unsigned& fld)
  {
#ifdef PARANOID
    if (fld != 0)
    {
      std::stringstream error_stream;
      error_stream << "HeleShaw elements only store a single field so fld "
                      "must be 0 rather"
                   << " than " << fld << std::endl;
      throw OomphLibError(error_stream.str(),
                          "ProjectableHeleShawElement::data_values_of_field()",
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Create the vector
    Vector<std::pair<Data*, unsigned>> data_values;

    // Loop over all nodes
    unsigned nnod = this->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      // Add the data value associated field: The node itself
      data_values.push_back(std::make_pair(this->node_pt(j), fld));
    }

    // Return the vector
    return data_values;
  }

  /// \short Number of fields to be projected: Just one
  template<class HELE_SHAW_ELEMENT>
  unsigned ProjectableHeleShawElement<
    HELE_SHAW_ELEMENT>::nfields_for_projection()
  {
    return 1;
  }

  /// \short Number of history values to be stored for fld-th field.
  template<class HELE_SHAW_ELEMENT>
  unsigned ProjectableHeleShawElement<
    HELE_SHAW_ELEMENT>::nhistory_values_for_projection(const unsigned& fld)
  {
#ifdef PARANOID
    if (fld != 0)
    {
      std::stringstream error_stream;
      error_stream << "HeleShaw elements only store a single field so fld "
                      "must be 0 rather"
                   << " than " << fld << std::endl;
      throw OomphLibError(
        error_stream.str(),
        "ProjectableHeleShawElement::nhistory_values_for_projection()",
        OOMPH_EXCEPTION_LOCATION);
    }
#endif
    return this->node_pt(0)->ntstorage();
  }

  ///\short Number of positional history values
  template<class HELE_SHAW_ELEMENT>
  unsigned ProjectableHeleShawElement<
    HELE_SHAW_ELEMENT>::nhistory_values_for_coordinate_projection()
  {
    return this->node_pt(0)->position_time_stepper_pt()->ntstorage();
  }

  /// \short Return Jacobian of mapping and shape functions of field fld
  /// at local coordinate s
  template<class HELE_SHAW_ELEMENT>
  double ProjectableHeleShawElement<
    HELE_SHAW_ELEMENT>::jacobian_and_shape_of_field(const unsigned& fld,
                                                    const Vector<double>& s,
                                                    Shape& psi)
  {
#ifdef PARANOID
    if (fld != 0)
    {
      std::stringstream error_stream;
      error_stream << "HeleShaw elements only store a single field so fld "
                      "must be 0 rather"
                   << " than " << fld << std::endl;
      throw OomphLibError(
        error_stream.str(),
        "ProjectableHeleShawElement::jacobian_and_shape_of_field()",
        OOMPH_EXCEPTION_LOCATION);
    }
#endif
    unsigned n_dim = this->dim();
    unsigned n_node = this->nnode();
    Shape test(n_node);
    DShape dpsidx(n_node, n_dim), dtestdx(n_node, n_dim);
    double J =
      this->dshape_and_dtest_eulerian_hele_shaw(s, psi, dpsidx, test, dtestdx);
    return J;
  }

  /// \short Return interpolated field fld at local coordinate s, at time
  /// level t (t=0: present; t>0: history values)
  template<class HELE_SHAW_ELEMENT>
  double ProjectableHeleShawElement<HELE_SHAW_ELEMENT>::get_field(
    const unsigned& t, const unsigned& fld, const Vector<double>& s)
  {
#ifdef PARANOID
    if (fld != 0)
    {
      std::stringstream error_stream;
      error_stream << "HeleShaw elements only store a single field so fld "
                      "must be 0 rather"
                   << " than " << fld << std::endl;
      throw OomphLibError(error_stream.str(),
                          "ProjectableHeleShawElement::jget_field()",
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    // Find the index at which the variable is stored
    unsigned p_nodal_index = this->p_index_hele_shaw();

    // Local shape function
    unsigned n_node = this->nnode();
    Shape psi(n_node);

    // Find values of shape function
    this->shape(s, psi);

    // Initialise value of u
    double interpolated_p = 0.0;

    // Sum over the local nodes
    for (unsigned l = 0; l < n_node; l++)
    {
      interpolated_p += this->nodal_value(l, p_nodal_index) * psi[l];
    }
    return interpolated_p;
  }

  /// Return number of values in field fld: One per node
  template<class HELE_SHAW_ELEMENT>
  unsigned ProjectableHeleShawElement<HELE_SHAW_ELEMENT>::nvalue_of_field(
    const unsigned& fld)
  {
#ifdef PARANOID
    if (fld != 0)
    {
      std::stringstream error_stream;
      error_stream << "HeleShaw elements only store a single field so fld "
                      "must be 0 rather"
                   << " than " << fld << std::endl;
      throw OomphLibError(error_stream.str(),
                          "ProjectableHeleShawElement::nvalue_of_field()",
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    return this->nnode();
  }

  /// Return local equation number of value j in field fld.
  template<class HELE_SHAW_ELEMENT>
  int ProjectableHeleShawElement<HELE_SHAW_ELEMENT>::local_equation(
    const unsigned& fld, const unsigned& j)
  {
#ifdef PARANOID
    if (fld != 0)
    {
      std::stringstream error_stream;
      error_stream << "HeleShaw elements only store a single field so fld "
                      "must be 0 rather"
                   << " than " << fld << std::endl;
      throw OomphLibError(error_stream.str(),
                          "ProjectableHeleShawElement::local_equation()",
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    const unsigned p_nodal_index = this->p_index_hele_shaw();
    return this->nodal_local_eqn(j, p_nodal_index);
  }

  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<ProjectableHeleShawElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };

  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as
  /// that for the underlying wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<ProjectableHeleShawElement<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };


} // namespace oomph
#endif
