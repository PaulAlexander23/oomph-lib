#ifndef my_element_header
#define my_element_header

#include "generic.h"
#include "axisym_navier_stokes.h"

namespace oomph
{
  class MyElement : public ElementWithError<AxisymmetricTTaylorHoodElement>
  {
  private:
    unsigned Region_id;

  public:
    void set_region_id(const unsigned& value)
    {
      Region_id = value;
    }

    const unsigned get_region_id()
    {
      return Region_id;
    }

    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates
      Vector<double> s(2);

      // Tecplot header info
      outfile << tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        get_s_plot(iplot, nplot, s);

        // Coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << interpolated_x(s, i) << " ";
        }

        // Velocities
        for (unsigned i = 0; i < 3; i++)
        {
          outfile << interpolated_u_axi_nst(s, i) << " ";
        }

        // Pressure
        outfile << interpolated_p_axi_nst(s) << " ";

        // Output the stored error
        outfile << this->get_error() << " ";

        // Output the element size
        outfile << this->size() << " ";

        // Output the continuity residual
        const double r = interpolated_x(s, 0);
        const double cont_res = this->interpolated_dudx_axi_nst(s, 0, 0) +
                                1 / r * this->interpolated_u_axi_nst(s, 0) +
                                this->interpolated_dudx_axi_nst(s, 1, 1);
        outfile << cont_res << " ";

        outfile << std::endl;
      }
      outfile << std::endl;

      // Write tecplot footer (e.g. FE connectivity lists)
      write_tecplot_zone_footer(outfile, nplot);
    }
  };


  //=======================================================================
  /// Face geometry of the 2D Taylor_Hood elements
  //=======================================================================
  template<>
  class FaceGeometry<MyElement> : public virtual TElement<1, 3>
  {
  public:
    /// Constructor: Call constructor of base
    FaceGeometry() : TElement<1, 3>() {}
  };
} // namespace oomph

#endif
