// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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

// Oomph-lib includes
#include "generic.h"
#include "hele_shaw.h"

using namespace std;
using namespace oomph;

//===start_of_main======================================================
/// Driver code: Testing Hele-Shaw parameter converter
//======================================================================
int main(int argc, char* argv[])
{
  std::ofstream out_stream;
  out_stream.open("OUTPUT");

  double length_ratio;
  double depth_ratio;
  double time_ratio;
  double pressure_ratio;
  get_ratios_of_q_nd_to_capillary_nd(
    length_ratio, depth_ratio, time_ratio, pressure_ratio);
  out_stream << length_ratio << ", ";
  out_stream << depth_ratio << ", ";
  out_stream << time_ratio << ", ";
  out_stream << pressure_ratio << ", ";

  double Q = 0.05;
  double V = 0.665;
  double h = 0.024;
  double Ca;
  double new_V;
  double new_h;
  convert_parameters_from_q_nd_to_capillary_nd(Q, V, h, Ca, new_V, new_h);
  out_stream << Ca << ", ";
  out_stream << new_V << ", ";
  out_stream << new_h << ", ";
  convert_parameters_from_capillary_nd_to_q_nd(Ca, new_V, new_h, Q, V, h);
  out_stream << Q << ", ";
  out_stream << V << ", ";
  out_stream << h << ", ";

  out_stream << endl;
  out_stream.close();

  return (EXIT_SUCCESS);
} // end_of_main
