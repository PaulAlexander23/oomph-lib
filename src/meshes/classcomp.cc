#include "classcomp.h"

#ifdef OOMPH_HAS_MPI

// ===================================================================
// The comparison class for the map that sorts the nodes on the
// shared boundary (using a lexicographic order)
// ===================================================================
// Comparison operator for "lower left" ordering
bool classcomp::operator()(const std::pair<double, double>& lhs,
                           const std::pair<double, double>& rhs) const
{
  double diff_y = lhs.second - rhs.second;
  if (diff_y < -Tol) // (lhs.second < rhs.second)
  {
    return true;
  }
  else
  {
    // Are they "equal" with 1.0e-14 tolerance?
    if (diff_y < Tol) // (lhs.second == rhs.second)
    {
#ifdef PARANOID
      double diff_x = lhs.first - rhs.first;
      if (fabs(diff_x) < Tol)
      {
        std::ostringstream warning_message;
        warning_message
          << "Dodgy \"lower left\" (lexicographic) comparison "
          << "of points with cooordinates: "
          << " lhs = ( " << lhs.first << " , " << lhs.second << " ) \n"
          << " rhs = ( " << rhs.first << " , " << rhs.second << " ) \n"
          << "x and y coordinates differ by less than tolerance!\n"
          << "diff_x = " << diff_x << "\n"
          << "diff_y = " << diff_y << "\n"
          << "Tol    = " << Tol << "\n";
        OomphLibError(warning_message.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
      }
#endif
      if (lhs.first < rhs.first)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }


  // if (lhs.second < rhs.second)
  //  {
  //   return true;
  //  }
  // else
  //  {
  //   // // Are "equal" with 1.0e-14 tolerance
  //   // if (lhs.second - rhs.second < 1.0e-14)
  //   // Are equal?
  //   if (lhs.second == rhs.second)
  //    {
  //     if (lhs.first < rhs.first)
  //      {
  //       return true;
  //      }
  //     else
  //      {
  //       return false;
  //      }
  //    }
  //   else
  //    {
  //     return false;
  //    }
  //  }
}

// Assign value for tolerance
double classcomp::Tol = 1.0e-14;

classcomp Bottom_left_sorter;

#endif
