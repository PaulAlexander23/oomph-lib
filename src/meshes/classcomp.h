#ifndef OOMPH_CLASSCOMP_HEADER
#define OOMPH_CLASSCOMP_HEADER

#include <utility>

#ifdef OOMPH_HAS_MPI

// ===================================================================
// The comparison class for the map that sorts the nodes on the
// shared boundary (using a lexicographic order)
// ===================================================================
struct classcomp
{
  // Tolerance for lower-left comparison
  static double Tol;


  // Comparison operator for "lower left" ordering
  bool operator()(const std::pair<double, double>& lhs,
                  const std::pair<double, double>& rhs) const;
}; // struct classcomp

#endif

#endif
