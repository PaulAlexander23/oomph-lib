#ifndef OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER
#define OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "integral.h"
#include "info_element.h"

#include "relaxing_bubble_parameters.h"
#include "my_constraint_elements.h"
#include "spatiotemporal_tolerences.h"

namespace oomph
{
  template<class ELEMENT>
  class RelaxingBubbleProblem : public Problem
  {
    /// Variables
  public:
    /// Enumeration of mesh boundaries
    enum
    {
      Left_boundary_id = 0,
      Top_boundary_id = 1,
      Right_boundary_id = 2,
      Bottom_boundary_id = 3,
      First_bubble_boundary_id = 4,
      Second_bubble_boundary_id = 5,
    };

  private:
    TriangleMeshPolygon* Outer_boundary_polygon_pt;
    TriangleMeshPolygon* Surface_polygon_pt;

    RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;
    Mesh* Surface_mesh_pt;
    Mesh* Volume_mesh_pt;
    Mesh* X_mom_mesh_pt;
    // Mesh* Y_mom_mesh_pt;
    Mesh* Volume_constraint_mesh_pt;

    /// Pointer to the "surface" mesh
    Mesh* Inlet_surface_mesh_pt;
    Mesh* Outlet_surface_mesh_pt;
    Mesh* Flux_mesh_pt;

    Data* Volume_data_pt;
    Data* X_mom_data_pt;
    // Data* Y_mom_data_pt;
    Data* Bubble_pressure_data_pt;

    SpatiotemporalTolerances* Tolerances_pt;

    /// Functions
  public:
    /// Constructor
    RelaxingBubbleProblem();

    RelaxingBubbleProblem(SpatiotemporalTolerances* tolerances_pt);

    /// Destructor
    ~RelaxingBubbleProblem() {}

    /// Use a newton solve to set the initial conditions
    void solve_for_initial_conditions(DocInfo& doc_info);

    /// Iterate forward in time
    void iterate_timestepper(const double& t_step,
                             const double& t_final,
                             DocInfo& doc_info);

    /// Doc the solution
    void doc_solution(DocInfo& doc_info);

  private:
    void create_mesh();
    void create_data();
    void create_outer_boundary_polygon();
    void create_surface_polygon();
    void create_fluid_mesh();
    void create_surface_mesh();
    void create_volume_mesh();
    void create_x_com_mesh();
    // void create_y_com_mesh();
    void create_volume_constraint_mesh();
    void create_inlet_surface_mesh();
    void create_outlet_surface_mesh();
    void create_flux_mesh();

    void set_variable_and_function_pointers();

    void set_boundary_conditions();
    void fill_in_bubble_boundary_map(map<unsigned, bool>& is_on_bubble_bound);

    void actions_before_adapt();
    void actions_after_adapt();

    void doc_fluid_mesh(string filename);
    void doc_surface_mesh(string filename);
    void doc_volume(string filename);
    void doc_x_com(string filename);
    void doc_bubble_pressure(string filename);
    void doc_boundary(string filename);

    void compute_error_estimate(double& max_err, double& min_err);

    void delete_surface_mesh();
    void delete_volume_mesh();
    void delete_x_com_mesh();
    // void delete_y_com_mesh();
    void delete_volume_constraint_mesh();
    void delete_inlet_surface_mesh();
    void delete_outlet_surface_mesh();
    void delete_flux_mesh();
    void delete_mesh_pt(Mesh* mesh_pt);
  };

  template<class ELEMENT>
  RelaxingBubbleProblem<ELEMENT>::RelaxingBubbleProblem()
  {
    cout << "Constructor" << endl;
    Tolerances_pt = new SpatiotemporalTolerances;

    bool adaptive_timestepping = Tolerances_pt->get_adaptive_timestepping();
    add_time_stepper_pt(new BDF<2>(adaptive_timestepping));

    create_data();

    create_outer_boundary_polygon();
    create_surface_polygon();

    create_mesh();

    set_variable_and_function_pointers();

    set_boundary_conditions();

    // Setup equation numbering scheme
    cout << "Number of equations: " << endl;
    cout << assign_eqn_numbers() << endl;
    cout << "Number of unknowns: " << endl;
    cout << ndof() << endl;
  }

  template<class ELEMENT>
  RelaxingBubbleProblem<ELEMENT>::RelaxingBubbleProblem(
    SpatiotemporalTolerances* tolerances_pt)
  {
    cout << "Constructor" << endl;
    Tolerances_pt = tolerances_pt;

    bool adaptive_timestepping = Tolerances_pt->get_adaptive_timestepping();
    add_time_stepper_pt(new BDF<2>(adaptive_timestepping));

    create_data();

    create_outer_boundary_polygon();
    create_surface_polygon();

    create_mesh();

    set_variable_and_function_pointers();

    set_boundary_conditions();

    // Setup equation numbering scheme
    cout << "Number of equations: " << endl;
    cout << assign_eqn_numbers() << endl;
    cout << "Number of unknowns: " << endl;
    cout << ndof() << endl;
  }

  /// Use a newton solve to set the initial conditions
  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::solve_for_initial_conditions(
    DocInfo& doc_info)
  {
    cout << "Solve for initial conditions" << endl;
    cout << "Steady Newton solve" << endl;
    // Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
    // const unsigned max_adapt = 1;
    // steady_newton_solve(max_adapt);

    double dt = Tolerances_pt->get_initial_timestep();
    initialise_dt(dt);
    assign_initial_values_impulsive(dt);
    Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

    doc_solution(doc_info);
  };

  /// Iterate forward in time
  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::iterate_timestepper(
    const double& t_step, const double& t_final, DocInfo& doc_info)
  {
    cout << "Iterate timestepper" << endl;


    unsigned n_timestep = ceil(t_final / t_step);

    unsigned max_adapt = 0;
    bool is_first_step = true;
    unsigned adapt_interval = Tolerances_pt->get_remesh_interval();
    for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++)
    {
      cout << "t: " << time() << endl;

      // max_newton_iterations() = 1;
      // newton_solver_tolerance() = 1e10;
      if (i_timestep % adapt_interval == 3 && !is_first_step)
      {
        max_adapt = 1;
      }
      else
      {
        max_adapt = 0;
      }

      // DoubleVector residuals;
      // DenseDoubleMatrix jacobian;
      // DoubleVector residualsFD;
      // DenseDoubleMatrix jacobianFD(ndof());

      // cout << "Get Jacobian" << endl;
      // get_jacobian(residuals, jacobian);
      ////get_fd_jacobian(residualsFD, jacobianFD);

      // for (unsigned i = 0; i < 1280; i++)
      //{
      //  printf("i: %4u", i);
      //  for (unsigned j = 1280 - 2; j < 1280; j++)
      //  {
      //    printf(", act: %8.5f, exp: %8.5f", jacobian(j, i), jacobianFD(j,
      //    i));
      //  }
      //  cout << endl;
      //}

      // double j = 0;
      // cout << "Residuals" << endl;
      // for (unsigned i = 0; i < ndof(); i++)
      // {
      //   j = residuals[i];
      //   cout << "i: " << i << ", j: " << j << endl;
      // }

      /// Allow (significant??) adapting on the first step
      if (is_first_step && Tolerances_pt->get_remesh_initial_condition())
      {
        max_adapt = 1;
      }

      cout << "Unsteady Newton solve" << endl;
      unsteady_newton_solve(t_step, max_adapt, is_first_step);
      Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

      /// Set the first step flag to false
      if (is_first_step)
      {
        is_first_step = false;
      }

      doc_solution(doc_info);
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_outer_boundary_polygon()
  {
    double domain_length = 1.0;
    double domain_width = 1.0;

    /// Create a mesh curve section for each boundary
    Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

    // Each polyline only has two vertices -- provide storage for their
    // coordinates
    Vector<Vector<double>> vertex_coord(2);
    for (unsigned i = 0; i < 2; i++)
    {
      vertex_coord[i].resize(2);
    }

    // First polyline: bottom
    vertex_coord[0][0] = 0.0;
    vertex_coord[0][1] = 0.0;
    vertex_coord[1][0] = domain_length;
    vertex_coord[1][1] = 0.0;

    // Build the 1st boundary polyline
    boundary_polyline_pt[0] =
      new TriangleMeshPolyLine(vertex_coord, Bottom_boundary_id);

    // Second boundary polyline: right
    vertex_coord[0][0] = vertex_coord[1][0];
    vertex_coord[0][1] = vertex_coord[1][1];
    vertex_coord[1][0] = domain_length;
    vertex_coord[1][1] = domain_width;

    // Build the 2nd boundary polyline
    boundary_polyline_pt[1] =
      new TriangleMeshPolyLine(vertex_coord, Right_boundary_id);

    // Third boundary polyline: top
    vertex_coord[0][0] = vertex_coord[1][0];
    vertex_coord[0][1] = vertex_coord[1][1];
    vertex_coord[1][0] = 0.0;
    vertex_coord[1][1] = domain_width;

    // Build the 3rd boundary polyline
    boundary_polyline_pt[2] =
      new TriangleMeshPolyLine(vertex_coord, Top_boundary_id);

    // Fourth boundary polyline: left
    vertex_coord[0][0] = vertex_coord[1][0];
    vertex_coord[0][1] = vertex_coord[1][1];
    vertex_coord[1][0] = 0.0;
    vertex_coord[1][1] = 0.0;

    // Build the 4th boundary polyline
    boundary_polyline_pt[3] =
      new TriangleMeshPolyLine(vertex_coord, Left_boundary_id);

    // Create the triangle mesh polygon for outer boundary
    Outer_boundary_polygon_pt = new TriangleMeshPolygon(boundary_polyline_pt);
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_surface_polygon()
  {
    double x_center = relaxing_bubble::bubble_initial_centre_x;
    double y_center = relaxing_bubble::bubble_initial_centre_y;
    double minor_radius =
      (relaxing_bubble::target_bubble_volume) /
      (MathematicalConstants::Pi * relaxing_bubble::major_radius);
    unsigned npoints =
      Tolerances_pt->get_initial_number_of_polynomial_vertices();
    double zeta_step = MathematicalConstants::Pi / double(npoints - 1);

    // Intrinsic coordinate along GeomObject defining the bubble
    Vector<double> zeta(1);

    // Position vector to GeomObject defining the bubble
    Vector<double> coord(2);

    Ellipse* ellipse_pt =
      new Ellipse(relaxing_bubble::major_radius, minor_radius);

    /// Triangle requires at least two polylines for a single polygon
    Vector<TriangleMeshCurveSection*> surface_polyline_pt(2);

    // Vertex coordinates
    Vector<Vector<double>> bubble_vertex(npoints);

    for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
    {
      bubble_vertex[ipoint].resize(2);

      // Get the coordinates
      zeta[0] = zeta_step * double(ipoint);
      ellipse_pt->position(zeta, coord);

      // Shift
      bubble_vertex[ipoint][0] = coord[0] + x_center;
      bubble_vertex[ipoint][1] = coord[1] + y_center;
    }

    // Build the 1st bubble polyline
    surface_polyline_pt[0] =
      new TriangleMeshPolyLine(bubble_vertex, First_bubble_boundary_id);

    for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
    {
      bubble_vertex[ipoint].resize(2);

      // Get the coordinates
      zeta[0] = zeta_step * double(ipoint) + MathematicalConstants::Pi;
      ellipse_pt->position(zeta, coord);

      // Shift
      bubble_vertex[ipoint][0] = coord[0] + x_center;
      bubble_vertex[ipoint][1] = coord[1] + y_center;
    }

    // Build the 2nd bubble polyline
    surface_polyline_pt[1] =
      new TriangleMeshPolyLine(bubble_vertex, Second_bubble_boundary_id);

    // Define coordinates of a point inside the bubble
    Vector<double> bubble_center(2);
    bubble_center[0] = x_center;
    bubble_center[1] = y_center;

    double max_length = Tolerances_pt->get_maximum_polyline_segment_length();
    surface_polyline_pt[0]->set_maximum_length(max_length);
    surface_polyline_pt[1]->set_maximum_length(max_length);

    // Create closed polygon from two polylines
    Surface_polygon_pt =
      new TriangleMeshPolygon(surface_polyline_pt, bubble_center);

    Surface_polygon_pt->set_polyline_refinement_tolerance(
      Tolerances_pt->get_polyline_refinement_tolerance());
    Surface_polygon_pt->set_polyline_unrefinement_tolerance(
      Tolerances_pt->get_polyline_unrefinement_tolerance());
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_data()
  {
    Volume_data_pt = new Data(1);
    X_mom_data_pt = new Data(1);
    // Y_mom_data_pt = new Data(1);

    Bubble_pressure_data_pt = new Data(1);
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_mesh()
  {
    create_fluid_mesh();
    Surface_mesh_pt = new Mesh;
    create_surface_mesh();
    Volume_mesh_pt = new Mesh;
    create_volume_mesh();
    X_mom_mesh_pt = new Mesh;
    create_x_com_mesh();
    // Y_mom_mesh_pt = new Mesh;
    // create_y_com_mesh();
    Volume_constraint_mesh_pt = new Mesh;
    create_volume_constraint_mesh();
    Flux_mesh_pt = new Mesh;
    create_flux_mesh();
    Inlet_surface_mesh_pt = new Mesh;
    create_inlet_surface_mesh();
    Outlet_surface_mesh_pt = new Mesh;
    create_outlet_surface_mesh();

    add_sub_mesh(Fluid_mesh_pt);
    add_sub_mesh(Surface_mesh_pt);
    add_sub_mesh(Volume_mesh_pt);
    add_sub_mesh(X_mom_mesh_pt);
    // add_sub_mesh(Y_mom_mesh_pt);
    add_sub_mesh(Volume_constraint_mesh_pt);
    add_sub_mesh(Inlet_surface_mesh_pt);
    add_sub_mesh(Outlet_surface_mesh_pt);
    add_sub_mesh(Flux_mesh_pt);

    build_global_mesh();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_fluid_mesh()
  {
    // Target area for initial mesh
    double uniform_element_area =
      Tolerances_pt->get_initial_target_element_area();

    TriangleMeshClosedCurve* outer_closed_curve_pt = Outer_boundary_polygon_pt;
    Vector<TriangleMeshClosedCurve*> surface_closed_curve_pt(1);
    surface_closed_curve_pt[0] = Surface_polygon_pt;

    // Use the TriangleMeshParameters object for gathering all
    // the necessary arguments for the TriangleMesh object
    TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

    // Define the holes on the boundary
    triangle_mesh_parameters.internal_closed_curve_pt() =
      surface_closed_curve_pt;

    // Define the maximum element areas
    triangle_mesh_parameters.element_area() = uniform_element_area;

    // Create the mesh
    Fluid_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
      triangle_mesh_parameters, this->time_stepper_pt());

    // Set error estimator for fluid mesh
    Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
    Fluid_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

    // Set targets for spatial adaptivity
    Fluid_mesh_pt->max_permitted_error() =
      Tolerances_pt->get_maximum_permitted_error();
    Fluid_mesh_pt->min_permitted_error() =
      Tolerances_pt->get_minimum_permitted_error();
    Fluid_mesh_pt->max_element_size() =
      Tolerances_pt->get_maximum_element_size();
    Fluid_mesh_pt->min_element_size() =
      Tolerances_pt->get_minimum_element_size();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_surface_mesh()
  {
    for (unsigned i_boundary = First_bubble_boundary_id;
         i_boundary < Second_bubble_boundary_id + 1;
         i_boundary++)
    {
      unsigned n_element = Fluid_mesh_pt->nboundary_element(i_boundary);
      for (unsigned i_element = 0; i_element < n_element; i_element++)
      {
        ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(i_boundary, i_element));

        // Find the index of the face of element e along boundary b
        int face_index =
          Fluid_mesh_pt->face_index_at_boundary(i_boundary, i_element);

        HeleShawInterfaceElement<ELEMENT>* interface_element_pt =
          new HeleShawInterfaceElement<ELEMENT>(fluid_element_pt, face_index);

        // Add the appropriate boundary number
        interface_element_pt->set_boundary_number_in_bulk_mesh(i_boundary);

        // Add it to the mesh
        Surface_mesh_pt->add_element_pt(interface_element_pt);
      }
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_volume_mesh()
  {
    InfoElement* volume_element_pt = new InfoElement;
    Volume_mesh_pt->add_element_pt(volume_element_pt);
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_x_com_mesh()
  {
    InfoElement* x_mom_element_pt = new InfoElement;
    X_mom_mesh_pt->add_element_pt(x_mom_element_pt);
  }

  // template<class ELEMENT>
  // void RelaxingBubbleProblem<ELEMENT>::create_y_com_mesh()
  //{
  //  InfoElement* y_mom_element_pt = new InfoElement;
  //  Y_mom_mesh_pt->add_element_pt(y_mom_element_pt);
  //}

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_volume_constraint_mesh()
  {
    MyConstraintElement* vol_constraint_element =
      new MyConstraintElement(&relaxing_bubble::target_fluid_volume,
                              Volume_data_pt->value_pt(0),
                              Bubble_pressure_data_pt,
                              0);
    Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_inlet_surface_mesh()
  {
    const unsigned boundary = Left_boundary_id;

    unsigned n_element = this->Fluid_mesh_pt->nboundary_element(boundary);
    for (unsigned n = 0; n < n_element; n++)
    {
      ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
        this->Fluid_mesh_pt->boundary_element_pt(boundary, n));

      int face_index = this->Fluid_mesh_pt->face_index_at_boundary(boundary, n);

      HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
        new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
          fluid_element_pt,
          face_index,
          this->Flux_mesh_pt->element_pt(0)->internal_data_pt(0));

      this->Inlet_surface_mesh_pt->add_element_pt(flux_element_pt);
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_outlet_surface_mesh()
  {
    const unsigned boundary = Right_boundary_id;

    unsigned n_element = this->Fluid_mesh_pt->nboundary_element(boundary);
    for (unsigned n = 0; n < n_element; n++)
    {
      ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
        this->Fluid_mesh_pt->boundary_element_pt(boundary, n));

      int face_index = this->Fluid_mesh_pt->face_index_at_boundary(boundary, n);

      HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
        new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
          fluid_element_pt,
          face_index,
          this->Flux_mesh_pt->element_pt(1)->internal_data_pt(0));

      this->Outlet_surface_mesh_pt->add_element_pt(flux_element_pt);
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_flux_mesh()
  {
    unsigned number_of_values = 1;
    Data* inlet_integral_data_pt = new Data(number_of_values);
    this->Flux_mesh_pt->add_element_pt(new InfoElement(inlet_integral_data_pt));
    Data* outlet_integral_data_pt = new Data(number_of_values);
    this->Flux_mesh_pt->add_element_pt(
      new InfoElement(outlet_integral_data_pt));
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::set_variable_and_function_pointers()
  {
    bool fd_jacobian = true;

    relaxing_bubble::bubble_pressure_pt = Bubble_pressure_data_pt->value_pt(0);

    /// Set fluid mesh function pointers
    unsigned n_element = Fluid_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast from GeneralisedElement to the present element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

      // Set the constitutive law for pseudo-elastic mesh deformation
      el_pt->constitutive_law_pt() = relaxing_bubble::constitutive_law_pt;
      el_pt->upper_wall_fct_pt() = relaxing_bubble::upper_wall_fct;
      el_pt->add_external_data(Volume_data_pt, fd_jacobian);
      // el_pt->add_volume_data_pt(Volume_data_pt);
      el_pt->add_external_data(X_mom_data_pt, fd_jacobian);
      // el_pt->add_external_data(Y_mom_data_pt, fd_jacobian);
      el_pt->add_external_data(Bubble_pressure_data_pt, fd_jacobian);
    }

    n_element = Surface_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      HeleShawInterfaceElement<ELEMENT>* interface_element_pt =
        dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
          Surface_mesh_pt->element_pt(e));
      interface_element_pt->ca_inv_pt() = &relaxing_bubble::ca_inv;
      interface_element_pt->st_pt() = &relaxing_bubble::st;
      interface_element_pt->aspect_ratio_pt() = &relaxing_bubble::alpha;
      interface_element_pt->upper_wall_fct_pt() = relaxing_bubble::upper_wall_fct;
      interface_element_pt->wall_speed_fct_pt() = relaxing_bubble::wall_speed_fct;
      interface_element_pt->bubble_pressure_fct_pt() =
        relaxing_bubble::bubble_pressure_fct;

      interface_element_pt->add_external_data(Volume_data_pt, fd_jacobian);
      interface_element_pt->add_external_data(Bubble_pressure_data_pt,
                                              fd_jacobian);
    }

    unsigned i_element = 0;
    InfoElement* volume_element_pt =
      dynamic_cast<InfoElement*>(Volume_mesh_pt->element_pt(i_element));
    volume_element_pt->add_data_pt(Volume_data_pt);

    InfoElement* x_mom_element_pt =
      dynamic_cast<InfoElement*>(X_mom_mesh_pt->element_pt(i_element));
    x_mom_element_pt->add_data_pt(X_mom_data_pt);

    // InfoElement* y_mom_element_pt =
    //  dynamic_cast<InfoElement*>(Y_mom_mesh_pt->element_pt(i_element));
    // y_mom_element_pt->add_data_pt(Y_mom_data_pt);

    MyConstraintElement* vol_constraint_element =
      dynamic_cast<MyConstraintElement*>(
        Volume_constraint_mesh_pt->element_pt(i_element));

    vol_constraint_element->add_external_data(Volume_data_pt, fd_jacobian);

    /// Info mesh
    unsigned index = 0;
    relaxing_bubble::inlet_b3_pt =
      this->Flux_mesh_pt->element_pt(0)->internal_data_pt(0)->value_pt(index);
    relaxing_bubble::outlet_b3_pt =
      this->Flux_mesh_pt->element_pt(1)->internal_data_pt(0)->value_pt(index);

    /// Inlet mesh
    // Find number of elements in mesh
    n_element = this->Inlet_surface_mesh_pt->nelement();

    // Loop over the elements to set up element-specific
    // things that cannot be handled by constructor
    for (unsigned i = 0; i < n_element; i++)
    {
      // Upcast from GeneralElement to the present element
      HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
        dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
          this->Inlet_surface_mesh_pt->element_pt(i));

      // Set the Neumann function pointer
      el_pt->flux_fct_pt() = &relaxing_bubble::get_inlet_flux_bc;
    }

    /// Outlet mesh
    // Find number of elements in mesh
    n_element = this->Outlet_surface_mesh_pt->nelement();

    // Loop over the elements to set up element-specific
    // things that cannot be handled by constructor
    for (unsigned i = 0; i < n_element; i++)
    {
      // Upcast from GeneralElement to the present element
      HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
        dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
          this->Outlet_surface_mesh_pt->element_pt(i));

      // Set the Neumann function pointer
      el_pt->flux_fct_pt() = &relaxing_bubble::get_outlet_flux_bc;
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::set_boundary_conditions()
  {
    // Map to record if a given boundary is on a bubble or not
    map<unsigned, bool> is_on_bubble_bound;
    fill_in_bubble_boundary_map(is_on_bubble_bound);


    const double fixed_pressure = 0.0;
    // Re-set the boundary conditions for fluid problem: All nodes are
    // free by default -- just pin the ones that have Dirichlet conditions
    // here.
    unsigned nbound = Fluid_mesh_pt->nboundary();
    for (unsigned ibound = 0; ibound < nbound; ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* node_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);

        // Pin pseudo-solid positions apart from bubble boundary which
        // we allow to move
        SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(node_pt);
        if (is_on_bubble_bound[ibound])
        {
          solid_node_pt->unpin_position(0);
          solid_node_pt->unpin_position(1);
        }
        else
        {
          solid_node_pt->pin_position(0);
          solid_node_pt->pin_position(1);

          node_pt->set_value(0, 0.0);
          node_pt->pin(0);
        }

        // if (ibound == Right_boundary_id)
        //{
        //  node_pt->set_value(0, fixed_pressure);
        //  node_pt->pin(0);
        //}
      }
    } // end loop over boundaries

    /// Pin tangential lagrange multiplier
    for (unsigned ibound = First_bubble_boundary_id;
         ibound < Second_bubble_boundary_id;
         ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* node_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
        ////                node_pt->pin(1); /// Normal lagrange multiplier
        ////                node_pt->pin(2); /// Curvature
        ////                node_pt->pin(3); /// Curvature
        node_pt->pin(4); /// Tangential lagrange multiplier
      }
    }

    Volume_data_pt->set_value(0, relaxing_bubble::target_fluid_volume);
    Volume_data_pt->unpin(0);

    X_mom_data_pt->set_value(0, 0.0);
    X_mom_data_pt->unpin(0);

    // Y_mom_data_pt->set_value(0, 0.0);
    // Y_mom_data_pt->unpin(0);

    // unsigned i_data = 0;
    // i_element = 0;
    // cout << "Get volume element" << endl;
    // InfoElement* volume_element_pt =
    //  dynamic_cast<InfoElement*>(Volume_mesh_pt->element_pt(i_element));
    // volume_element_pt->internal_data_pt(i_data)->set_value(0, 2.5);
    // volume_element_pt->internal_data_pt(i_data)->unpin(0);

    Bubble_pressure_data_pt->set_value(0, 16.0 / 3.0);
    Bubble_pressure_data_pt->unpin(0);

    /// Integral info mesh
    const unsigned value_index = 0;
    for (unsigned index = 0; index < 2; index++)
    {
      /// Integral value must be initialised to something non-zero?
      this->Flux_mesh_pt->element_pt(index)->internal_data_pt(0)->set_value(
        value_index, 1.0);
      this->Flux_mesh_pt->element_pt(index)->internal_data_pt(0)->unpin(
        value_index);
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::fill_in_bubble_boundary_map(
    map<unsigned, bool>& is_on_bubble_bound)
  {
    Vector<unsigned> bubble_bound_id =
      this->Surface_polygon_pt->polygon_boundary_id();
    // Get the number of boundary
    unsigned nbound = bubble_bound_id.size();
    // Fill in the map
    for (unsigned ibound = 0; ibound < nbound; ibound++)
    {
      // This boundary...
      unsigned bound_id = bubble_bound_id[ibound];
      // ...is on the bubble
      is_on_bubble_bound[bound_id] = true;
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
  {
    string doc_directory = doc_info.directory();

    doc_fluid_mesh(doc_directory + "soln" + to_string(doc_info.number()) +
                   ".dat");

    doc_surface_mesh(doc_directory + "surface" + to_string(doc_info.number()) +
                     ".dat");

    doc_volume(doc_directory + "volume.dat");

    doc_x_com(doc_directory + "com.dat");

    doc_bubble_pressure(doc_directory + "bubble_pressure.dat");

    // doc_boundary(doc_directory + "boundary.dat");

    doc_info.number()++;
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::actions_before_adapt()
  {
    cout << "Actions before adapt" << endl;

    delete_surface_mesh();
    delete_volume_mesh();
    delete_x_com_mesh();
    // delete_y_com_mesh();
    delete_volume_constraint_mesh();
    delete_flux_mesh();
    delete_inlet_surface_mesh();
    delete_outlet_surface_mesh();

    rebuild_global_mesh();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::actions_after_adapt()
  {
    cout << "Actions after adapt" << endl;

    Volume_data_pt = new Data(1);
    Bubble_pressure_data_pt = new Data(1);
    X_mom_data_pt = new Data(1);
    // Y_mom_data_pt = new Data(1);

    create_surface_mesh();
    create_volume_mesh();
    create_x_com_mesh();
    // create_y_com_mesh();
    create_volume_constraint_mesh();
    create_flux_mesh();
    create_inlet_surface_mesh();
    create_outlet_surface_mesh();

    rebuild_global_mesh();

    cout << "set_variable_and_function_pointers" << endl;
    set_variable_and_function_pointers();

    cout << "set_boundary_conditions" << endl;
    set_boundary_conditions();

    // Setup equation numbering scheme
    cout << "Number of equations: " << endl;
    cout << assign_eqn_numbers() << endl;
    cout << "Number of unknowns: " << endl;
    cout << ndof() << endl;
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::doc_fluid_mesh(string filename)
  {
    double max_err;
    double min_err;
    compute_error_estimate(max_err, min_err);

    ofstream output_stream;
    output_stream.open(filename);
    unsigned npoints = 3;
    Fluid_mesh_pt->output(output_stream, npoints);
    output_stream.close();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::doc_surface_mesh(string filename)
  {
    ofstream output_stream;
    output_stream.open(filename);

    unsigned n_element = Surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++)
    {
      dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
        Surface_mesh_pt->element_pt(e))
        ->output(output_stream, 3);
    }

    output_stream.close();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::doc_volume(string filename)
  {
    ofstream output_stream;
    output_stream.open(filename, ofstream::app);

    output_stream << time() << ", ";
    output_stream << Volume_data_pt->value(0) << endl;

    output_stream.close();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::doc_x_com(string filename)
  {
    ofstream output_stream;
    output_stream.open(filename, ofstream::app);

    output_stream << time() << ", ";
    output_stream << X_mom_data_pt->value(0) / Volume_data_pt->value(0) << endl;

    output_stream.close();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::doc_bubble_pressure(string filename)
  {
    ofstream output_stream;
    output_stream.open(filename, ofstream::app);

    output_stream << time() << ", ";
    output_stream << Bubble_pressure_data_pt->value(0) << endl;

    output_stream.close();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::doc_boundary(string filename)
  {
    // ofstream output_stream;
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::compute_error_estimate(double& max_err,
                                                              double& min_err)
  {
    // Get error estimator
    ErrorEstimator* err_est_pt = Fluid_mesh_pt->spatial_error_estimator_pt();

    // Get/output error estimates
    unsigned n_elements = Fluid_mesh_pt->nelement();
    Vector<double> elemental_error(n_elements);

    // We need a dynamic cast, get_element_errors from the Fluid_mesh_pt
    // Dynamic cast is used because get_element_errors require a Mesh* ans
    // not a SolidMesh*
    Mesh* fluid_mesh_pt = dynamic_cast<Mesh*>(Fluid_mesh_pt);
    err_est_pt->get_element_errors(fluid_mesh_pt, elemental_error);

    // Set errors for post-processing and find extrema
    max_err = 0.0;
    min_err = 1e6;
    for (unsigned e = 0; e < n_elements; e++)
    {
      dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e))
        ->set_error(elemental_error[e]);

      max_err = std::max(max_err, elemental_error[e]);
      min_err = std::min(min_err, elemental_error[e]);
    }

    std::cout << "Max error is " << max_err << std::endl;
    std::cout << "Min error is " << min_err << std::endl;
  }


  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_surface_mesh()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = Surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete surface element
      delete Surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Surface_mesh_pt->flush_element_and_node_storage();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_volume_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = Volume_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete Volume_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Volume_mesh_pt->flush_element_and_node_storage();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_x_com_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = X_mom_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete X_mom_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    X_mom_mesh_pt->flush_element_and_node_storage();
  }

  // template<class ELEMENT>
  // void RelaxingBubbleProblem<ELEMENT>::delete_y_com_mesh()
  //{
  //  // How many elements are in the volume mesh
  //  unsigned n_element = Y_mom_mesh_pt->nelement();

  //  // Loop over the elements
  //  for (unsigned e = 0; e < n_element; e++)
  //  {
  //    // Delete element
  //    delete Y_mom_mesh_pt->element_pt(e);
  //  }

  //  // Wipe the mesh
  //  Y_mom_mesh_pt->flush_element_and_node_storage();
  //}

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_volume_constraint_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = Volume_constraint_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete Volume_constraint_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Volume_constraint_mesh_pt->flush_element_and_node_storage();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_inlet_surface_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = Inlet_surface_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete Inlet_surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Inlet_surface_mesh_pt->flush_element_and_node_storage();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_outlet_surface_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = Outlet_surface_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete Outlet_surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Outlet_surface_mesh_pt->flush_element_and_node_storage();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_flux_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = Flux_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete Flux_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Flux_mesh_pt->flush_element_and_node_storage();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::delete_mesh_pt(Mesh* mesh_pt)
  {
    // How many elements are in the volume mesh
    unsigned n_element = mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    mesh_pt->flush_element_and_node_storage();
  }

} // namespace oomph
#endif
