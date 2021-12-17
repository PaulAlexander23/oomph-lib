#ifndef OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER
#define OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "integral.h"
#include "info_element.h"

#include "finger_parameters.h"
#include "my_constraint_elements.h"

namespace oomph
{
  template<class ELEMENT>
  class FingerProblem : public Problem
  {
    /// Variables
  public:
    /// Enumeration of mesh boundaries
    enum
    {
      Bottom_boundary_id = 0,
      Right_boundary_id = 1,
      Top_boundary_id = 2,
      Upper_inlet_boundary_id = 3,
      Finger_boundary_id = 4,
      Lower_inlet_boundary_id = 5,
    };

  private:
    TriangleMeshPolygon* Outer_boundary_polygon_pt;

    /// Fluid
    RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;

    Data* Volume_data_pt;
    Mesh* Volume_mesh_pt;

    Data* X_moment_data_pt;
    Mesh* X_moment_mesh_pt;

    /// Flux boundary
    Data* Top_inlet_integral_data_pt;
    Mesh* Top_inlet_surface_mesh_pt;

    Data* Bottom_inlet_integral_data_pt;
    Mesh* Bottom_inlet_surface_mesh_pt;

    Data* Outlet_integral_data_pt;
    Mesh* Outlet_surface_mesh_pt;

    /// Integral info from flux boundaries
    Mesh* Flux_mesh_pt;

    /// Finger
    // Mesh* Surface_mesh_pt;

    /// Fixed volume -> finger pressure
    // Data* Bubble_pressure_data_pt;
    // Mesh* Volume_constraint_mesh_pt;

    /// Functions
  public:
    /// Constructor
    FingerProblem();

    /// Destructor
    ~FingerProblem() {}

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
    void create_fluid_mesh();
    // void create_surface_mesh();
    void create_volume_mesh();
    void create_x_moment_mesh();
    // void create_volume_constraint_mesh();
    void create_flux_surface_meshes();
    void create_flux_mesh();

    void setup_variable_and_function_pointers();
    // void setup_surface_mesh_functions_pointers();

    void set_boundary_conditions();

    void actions_before_adapt();
    void actions_after_adapt();

    void doc_fluid_mesh(string filename);
    // void doc_surface_mesh(string filename);
    void doc_volume(string filename);
    // void doc_bubble_pressure(string filename);
    void doc_integrals(string filename);

    void compute_error_estimate(double& max_err, double& min_err);

    // void delete_surface_mesh();
    void delete_volume_mesh();
    void delete_x_moment_mesh();
    // void delete_volume_constraint_mesh();
    void delete_inlet_surface_mesh();
    void delete_flux_surface_meshes();
    void delete_flux_mesh();
    void delete_mesh_pt(Mesh* mesh_pt);
  };

  template<class ELEMENT>
  FingerProblem<ELEMENT>::FingerProblem()
  {
    cout << "Constructor" << endl;
    bool adaptive_timestepping = false;
    add_time_stepper_pt(new BDF<2>(adaptive_timestepping));

    create_data();

    create_outer_boundary_polygon();

    create_mesh();

    setup_variable_and_function_pointers();

    set_boundary_conditions();

    // Setup equation numbering scheme
    cout << "Number of equations: " << endl;
    cout << assign_eqn_numbers() << endl;
    cout << "Number of unknowns: " << endl;
    cout << ndof() << endl;
  }

  /// Use a newton solve to set the initial conditions
  template<class ELEMENT>
  void FingerProblem<ELEMENT>::solve_for_initial_conditions(DocInfo& doc_info)
  {
    cout << "Solve for initial conditions" << endl;
    cout << "Steady Newton solve" << endl;
    // Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
    // const unsigned max_adapt = 1;
    // steady_newton_solve(max_adapt);

    double dt = 5e-3;
    initialise_dt(dt);
    assign_initial_values_impulsive(dt);
    Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

    doc_solution(doc_info);
  };

  /// Iterate forward in time
  template<class ELEMENT>
  void FingerProblem<ELEMENT>::iterate_timestepper(const double& t_step,
                                                   const double& t_final,
                                                   DocInfo& doc_info)
  {
    cout << "Iterate timestepper" << endl;


    unsigned n_timestep = ceil(t_final / t_step);

    unsigned max_adapt = 0;
    bool is_first_step = true;
    unsigned adapt_interval = 5;
    for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++)
    {
      cout << "t: " << time() << endl;

      // max_newton_iterations() = 1;
      // newton_solver_tolerance() = 1e10;
      if (i_timestep % adapt_interval == adapt_interval - 1)
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
      // DenseDoubleMatrix jacobian_error(ndof());

      // cout << "Get Jacobian" << endl;
      // get_jacobian(residuals, jacobian);
      // cout << "Get FD Jacobian" << endl;
      // get_fd_jacobian(residualsFD, jacobianFD);

      // for (unsigned i = 0; i < ndof(); i++)
      //{
      //  for (unsigned j = 0; j < ndof(); j++)
      //  {
      //    jacobian_error(i, j) = jacobian(i, j) - jacobianFD(i, j);

      //    if (abs(jacobian_error(i, j)) > 1e-3)
      //    {
      //      printf("i: %4u, j: %4u, act: %8.5f, exp: %8.5f \n",
      //             i,
      //             j,
      //             jacobian(i, j),
      //             jacobianFD(i, j));
      //    }
      //  }
      //}

      // for (unsigned i = 0; i < ndof() + 1; i++)
      //{
      //  printf("i: %4u", i);
      //  for (unsigned j = ndof() - 4; j < ndof() + 1; j++)
      //  {
      //    // printf(", act: %9.2e", jacobian(j, i));
      //    printf(", act: %9.2e, exp: %9.2e", jacobian(j, i), jacobianFD(j,
      //    i));
      //  }
      //  cout << endl;
      //}

      // double j = 0;
      // cout << "Residuals" << endl;
      // for (unsigned i = 0; i < ndof() + 1; i++)
      //{
      //  j = residuals[i];
      //  cout << "i: " << i << ", j: " << j << endl;
      //}


      cout << "Unsteady Newton solve" << endl;
      unsteady_newton_solve(t_step, max_adapt, is_first_step);
      if (is_first_step)
      {
        is_first_step = false;
      }
      Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

      doc_solution(doc_info);
    }
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_outer_boundary_polygon()
  {
    double domain_length = 1.0;
    double domain_width = 1.0;

    /// Create a mesh curve section for each boundary
    Vector<TriangleMeshCurveSection*> boundary_polyline_pt(6);

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

    // Fourth boundary polyline: Upper inlet
    vertex_coord[0][0] = vertex_coord[1][0];
    vertex_coord[0][1] = vertex_coord[1][1];
    vertex_coord[1][0] = 0.0;
    vertex_coord[1][1] = domain_width * 0.5 + finger::finger_width * 0.5;

    // Build the 4th boundary polyline
    boundary_polyline_pt[3] =
      new TriangleMeshPolyLine(vertex_coord, Upper_inlet_boundary_id);


    // The finger polyline only has many vertices -- provide storage for their
    // coordinates
    double x_center = 0.0;
    double y_center = domain_width * 0.5;
    double major_radius = finger::finger_width / 2;
    double minor_radius = finger::finger_width / 2;
    Ellipse* ellipse_pt = new Ellipse(major_radius, minor_radius);
    unsigned n_point = 64;
    double zeta_step = MathematicalConstants::Pi / double(n_point - 1);
    // Intrinsic coordinate along GeomObject defining the bubble
    Vector<double> zeta(1);
    // Position vector to GeomObject defining the bubble
    Vector<double> coord(2);
    Vector<Vector<double>> finger_vertex_coord(n_point);
    for (unsigned i = 0; i < n_point; i++)
    {
      finger_vertex_coord[i].resize(2);
    }
    double y_step = finger::finger_width / double(n_point - 1);
    for (unsigned i_point = 0; i_point < n_point; i_point++)
    {
      // Get the coordinates
      zeta[0] = MathematicalConstants::Pi / 2 - zeta_step * double(i_point);
      ellipse_pt->position(zeta, coord);

      // Shift
      finger_vertex_coord[i_point][0] = coord[0] + x_center;
      finger_vertex_coord[i_point][1] = coord[1] + y_center;
    }

    // Build the 5th boundary polyline
    boundary_polyline_pt[4] =
      new TriangleMeshPolyLine(finger_vertex_coord, Finger_boundary_id);

    boundary_polyline_pt[4]->set_maximum_length(0.02);

    // Fourth boundary polyline: Lower inlet
    vertex_coord[0][0] = finger_vertex_coord[n_point - 1][0];
    vertex_coord[0][1] = finger_vertex_coord[n_point - 1][1];
    vertex_coord[1][0] = 0.0;
    vertex_coord[1][1] = 0.0;

    // Build the 6th boundary polyline
    boundary_polyline_pt[5] =
      new TriangleMeshPolyLine(vertex_coord, Lower_inlet_boundary_id);


    // Create the triangle mesh polygon for outer boundary
    Outer_boundary_polygon_pt = new TriangleMeshPolygon(boundary_polyline_pt);
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_data()
  {
    Volume_data_pt = new Data(1);
    X_moment_data_pt = new Data(1);

    // Bubble_pressure_data_pt = new Data(1);

    Outlet_integral_data_pt = new Data(1);
    Top_inlet_integral_data_pt = new Data(1);
    Bottom_inlet_integral_data_pt = new Data(1);
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_mesh()
  {
    create_fluid_mesh();
    // Surface_mesh_pt = new Mesh;
    // create_surface_mesh();
    Volume_mesh_pt = new Mesh;
    create_volume_mesh();
    X_moment_mesh_pt = new Mesh;
    create_x_moment_mesh();
    // Volume_constraint_mesh_pt = new Mesh;
    // create_volume_constraint_mesh();
    Flux_mesh_pt = new Mesh;
    create_flux_mesh();
    Outlet_surface_mesh_pt = new Mesh;
    Top_inlet_surface_mesh_pt = new Mesh;
    Bottom_inlet_surface_mesh_pt = new Mesh;
    create_flux_surface_meshes();

    add_sub_mesh(Fluid_mesh_pt);
    // add_sub_mesh(Surface_mesh_pt);
    add_sub_mesh(Volume_mesh_pt);
    add_sub_mesh(X_moment_mesh_pt);
    // add_sub_mesh(Volume_constraint_mesh_pt);
    add_sub_mesh(Outlet_surface_mesh_pt);
    add_sub_mesh(Top_inlet_surface_mesh_pt);
    add_sub_mesh(Bottom_inlet_surface_mesh_pt);
    add_sub_mesh(Flux_mesh_pt);

    build_global_mesh();
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_fluid_mesh()
  {
    // Target area for initial mesh
    double uniform_element_area = 1e-2;

    TriangleMeshClosedCurve* outer_closed_curve_pt = Outer_boundary_polygon_pt;

    // Use the TriangleMeshParameters object for gathering all
    // the necessary arguments for the TriangleMesh object
    TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

    // Define the maximum element areas
    triangle_mesh_parameters.element_area() = uniform_element_area;

    // Create the mesh
    Fluid_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
      triangle_mesh_parameters, this->time_stepper_pt());

    // Set error estimator for fluid mesh
    Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
    Fluid_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

    // Set targets for spatial adaptivity
    Fluid_mesh_pt->max_permitted_error() = 1e-3;
    Fluid_mesh_pt->min_permitted_error() = 5e-6;
    Fluid_mesh_pt->max_element_size() = 5e-2;
    Fluid_mesh_pt->min_element_size() = 1e-4;
  }

  // template<class ELEMENT>
  // void FingerProblem<ELEMENT>::create_surface_mesh()
  //{
  //  unsigned n_element = Fluid_mesh_pt->nboundary_element(Finger_boundary_id);
  //  cout << "Finger boundary elements: " << n_element << endl;
  //  for (unsigned i_element = 0; i_element < n_element; i_element++)
  //  {
  //    ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
  //      Fluid_mesh_pt->boundary_element_pt(Finger_boundary_id, i_element));

  //    // Find the index of the face of element e along boundary b
  //    int face_index =
  //      Fluid_mesh_pt->face_index_at_boundary(Finger_boundary_id, i_element);

  //    HeleShawInterfaceElement<ELEMENT>* interface_element_pt =
  //      new HeleShawInterfaceElement<ELEMENT>(fluid_element_pt, face_index);

  //    // Add the appropriate boundary number
  //    interface_element_pt->set_boundary_number_in_bulk_mesh(
  //      Finger_boundary_id);

  //    // Add it to the mesh
  //    Surface_mesh_pt->add_element_pt(interface_element_pt);
  //  }
  //}

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_volume_mesh()
  {
    InfoElement* volume_element_pt = new InfoElement;
    Volume_mesh_pt->add_element_pt(volume_element_pt);
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_x_moment_mesh()
  {
    InfoElement* x_moment_element_pt = new InfoElement;
    X_moment_mesh_pt->add_element_pt(x_moment_element_pt);
  }

  // template<class ELEMENT>
  // void FingerProblem<ELEMENT>::create_volume_constraint_mesh()
  //{
  //  MyConstraintElement* vol_constraint_element =
  //    new MyConstraintElement(finger::target_fluid_volume_fct,
  //                            Volume_data_pt->value_pt(0),
  //                            Bubble_pressure_data_pt,
  //                            0);
  //  Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);
  //}

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_flux_surface_meshes()
  {
    unsigned n_element =
      this->Fluid_mesh_pt->nboundary_element(Right_boundary_id);
    cout << "n_element:" << n_element << endl;
    for (unsigned n = 0; n < n_element; n++)
    {
      ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
        this->Fluid_mesh_pt->boundary_element_pt(Right_boundary_id, n));

      int face_index =
        this->Fluid_mesh_pt->face_index_at_boundary(Right_boundary_id, n);

      HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
        new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
          fluid_element_pt, face_index, Outlet_integral_data_pt);

      this->Outlet_surface_mesh_pt->add_element_pt(flux_element_pt);
    }

    n_element = this->Fluid_mesh_pt->nboundary_element(Upper_inlet_boundary_id);
    cout << "Upper boundary elements: " << n_element << endl;
    for (unsigned n = 0; n < n_element; n++)
    {
      ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
        this->Fluid_mesh_pt->boundary_element_pt(Upper_inlet_boundary_id, n));

      int face_index =
        this->Fluid_mesh_pt->face_index_at_boundary(Upper_inlet_boundary_id, n);

      HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
        new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
          fluid_element_pt, face_index, Top_inlet_integral_data_pt);

      this->Top_inlet_surface_mesh_pt->add_element_pt(flux_element_pt);
    }

    n_element = this->Fluid_mesh_pt->nboundary_element(Lower_inlet_boundary_id);
    cout << "Lower boundary elements: " << n_element << endl;
    for (unsigned n = 0; n < n_element; n++)
    {
      ELEMENT* fluid_element_pt = dynamic_cast<ELEMENT*>(
        this->Fluid_mesh_pt->boundary_element_pt(Lower_inlet_boundary_id, n));

      int face_index =
        this->Fluid_mesh_pt->face_index_at_boundary(Lower_inlet_boundary_id, n);

      HeleShawFluxElementWithInflowIntegral<ELEMENT>* flux_element_pt =
        new HeleShawFluxElementWithInflowIntegral<ELEMENT>(
          fluid_element_pt, face_index, Bottom_inlet_integral_data_pt);

      this->Bottom_inlet_surface_mesh_pt->add_element_pt(flux_element_pt);
    }
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::create_flux_mesh()
  {
    this->Flux_mesh_pt->add_element_pt(
      new InfoElement(Outlet_integral_data_pt));
    this->Flux_mesh_pt->add_element_pt(
      new InfoElement(Top_inlet_integral_data_pt));
    this->Flux_mesh_pt->add_element_pt(
      new InfoElement(Bottom_inlet_integral_data_pt));
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::setup_variable_and_function_pointers()
  {
    bool fd_jacobian = true;

    // finger::bubble_pressure_pt = Bubble_pressure_data_pt->value_pt(0);

    /// Set fluid mesh function pointers
    unsigned n_element = Fluid_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast from GeneralisedElement to the present element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

      // Set the constitutive law for pseudo-elastic mesh deformation
      el_pt->constitutive_law_pt() = finger::constitutive_law_pt;
      el_pt->upper_wall_fct_pt() = finger::upper_wall_fct;
      el_pt->add_external_data(Volume_data_pt);
      el_pt->add_external_data(X_moment_data_pt);
      // el_pt->add_external_data(Bubble_pressure_data_pt, fd_jacobian);
    }

    // setup_surface_mesh_functions_pointers();

    unsigned i_element = 0;
    InfoElement* volume_element_pt =
      dynamic_cast<InfoElement*>(Volume_mesh_pt->element_pt(i_element));
    volume_element_pt->add_data_pt(Volume_data_pt);

    InfoElement* x_moment_element_pt =
      dynamic_cast<InfoElement*>(X_moment_mesh_pt->element_pt(i_element));
    x_moment_element_pt->add_data_pt(X_moment_data_pt);

    //  MyConstraintElement* vol_constraint_element =
    //    dynamic_cast<MyConstraintElement*>(
    //      Volume_constraint_mesh_pt->element_pt(i_element));

    //  vol_constraint_element->add_external_data(Volume_data_pt, fd_jacobian);

    /// Info mesh
    unsigned index = 0;
    finger::outlet_b3_pt = Outlet_integral_data_pt->value_pt(0);
    finger::top_inlet_b3_pt = Top_inlet_integral_data_pt->value_pt(0);
    finger::bottom_inlet_b3_pt = Bottom_inlet_integral_data_pt->value_pt(0);


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
      el_pt->flux_fct_pt() = &finger::get_outlet_flux_bc;
      el_pt->add_external_data(Outlet_integral_data_pt, fd_jacobian);
    }

    // Find number of elements in mesh
    n_element = this->Top_inlet_surface_mesh_pt->nelement();

    // Loop over the elements to set up element-specific
    // things that cannot be handled by constructor
    for (unsigned i = 0; i < n_element; i++)
    {
      // Upcast from GeneralElement to the present element
      HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
        dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
          this->Top_inlet_surface_mesh_pt->element_pt(i));

      // Set the Neumann function pointer
      el_pt->flux_fct_pt() = &finger::get_top_inlet_flux_bc;
      el_pt->add_external_data(Top_inlet_integral_data_pt, fd_jacobian);
    }

    // Find number of elements in mesh
    n_element = this->Bottom_inlet_surface_mesh_pt->nelement();

    // Loop over the elements to set up element-specific
    // things that cannot be handled by constructor
    for (unsigned i = 0; i < n_element; i++)
    {
      // Upcast from GeneralElement to the present element
      HeleShawFluxElementWithInflowIntegral<ELEMENT>* el_pt =
        dynamic_cast<HeleShawFluxElementWithInflowIntegral<ELEMENT>*>(
          this->Bottom_inlet_surface_mesh_pt->element_pt(i));

      // Set the Neumann function pointer
      el_pt->flux_fct_pt() = &finger::get_bottom_inlet_flux_bc;
      el_pt->add_external_data(Bottom_inlet_integral_data_pt, fd_jacobian);
    }
  }

  // template<class ELEMENT>
  // void FingerProblem<ELEMENT>::setup_surface_mesh_functions_pointers()
  //{
  //  bool fd_jacobian = true;

  //  unsigned n_element = Surface_mesh_pt->nelement();
  //  for (unsigned e = 0; e < n_element; e++)
  //  {
  //    HeleShawInterfaceElement<ELEMENT>* interface_element_pt =
  //      dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
  //        Surface_mesh_pt->element_pt(e));
  //    interface_element_pt->ca_inv_pt() = &finger::ca_inv;
  //    interface_element_pt->st_pt() = &finger::st;
  //    interface_element_pt->aspect_ratio_pt() = &finger::alpha;
  //    interface_element_pt->upper_wall_fct_pt() = finger::upper_wall_fct;
  //    interface_element_pt->wall_speed_fct_pt() = finger::wall_speed_fct;
  //    interface_element_pt->bubble_pressure_fct_pt() =
  //      finger::bubble_pressure_fct;

  //    interface_element_pt->add_external_data(Volume_data_pt, fd_jacobian);
  //    interface_element_pt->add_external_data(Bubble_pressure_data_pt,
  //                                            fd_jacobian);
  //  }
  //}

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::set_boundary_conditions()
  {
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
        if (ibound == Finger_boundary_id)
        {
          // solid_node_pt->unpin_position(0);
          // solid_node_pt->unpin_position(1);
          solid_node_pt->pin_position(0);
          solid_node_pt->pin_position(1);
        }
        else
        {
          solid_node_pt->pin_position(0);
          solid_node_pt->pin_position(1);
        }

        // if (ibound == Upper_inlet_boundary_id ||
        //     ibound == Lower_inlet_boundary_id)
        // {
        //   node_pt->set_value(0, fixed_pressure);
        //   node_pt->pin(0);
        //   cout<<"pinning"<<endl;
        // }
      }
    } // end loop over boundaries

    Node* node_pt = Fluid_mesh_pt->boundary_node_pt(Right_boundary_id, 0);
    node_pt->set_value(0, fixed_pressure);
    node_pt->pin(0);

    /// Pin tangential lagrange multiplier
    // unsigned ibound = Finger_boundary_id;
    // unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
    // for (unsigned inod = 0; inod < num_nod; inod++)
    //{
    //  // Get node
    //  Node* node_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
    //  ////                node_pt->pin(1); /// Normal lagrange multiplier
    //  ////                node_pt->pin(2); /// Curvature
    //  ////                node_pt->pin(3); /// Curvature
    //  node_pt->pin(4); /// Tangential lagrange multiplier
    //}

    Volume_data_pt->set_value(0, finger::target_fluid_volume);
    Volume_data_pt->unpin(0);

    X_moment_data_pt->set_value(0, 0.0);
    X_moment_data_pt->unpin(0);

    // Bubble_pressure_data_pt->set_value(0, 0.0);
    // Bubble_pressure_data_pt->pin(0);

    Outlet_integral_data_pt->set_value(0, 1.0);
    Outlet_integral_data_pt->unpin(0);

    Top_inlet_integral_data_pt->set_value(0, 0.2);
    Top_inlet_integral_data_pt->unpin(0);

    Bottom_inlet_integral_data_pt->set_value(0, 0.2);
    Bottom_inlet_integral_data_pt->unpin(0);
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
  {
    string doc_directory = doc_info.directory();

    doc_fluid_mesh(doc_directory + "soln" + to_string(doc_info.number()) +
                   ".dat");
    // doc_surface_mesh(doc_directory + "surface" + doc_info.number() + ".dat");

    doc_volume(doc_directory + "volume.dat");

    // doc_bubble_pressure(doc_directory + "bubble_pressure.dat");

    doc_integrals(doc_directory + "integral.dat");

    doc_info.number()++;
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::actions_before_adapt()
  {
    cout << "Actions before adapt" << endl;

    // delete_surface_mesh();
    delete_volume_mesh();
    delete_x_moment_mesh();
    // delete_volume_constraint_mesh();
    delete_flux_surface_meshes();
    delete_flux_mesh();

    rebuild_global_mesh();
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::actions_after_adapt()
  {
    cout << "Actions after adapt" << endl;

    Volume_data_pt = new Data(1);
    X_moment_data_pt = new Data(1);

    // Bubble_pressure_data_pt = new Data(1);

    Outlet_integral_data_pt = new Data(1);
    Top_inlet_integral_data_pt = new Data(1);
    Bottom_inlet_integral_data_pt = new Data(1);

    // create_surface_mesh();
    create_volume_mesh();
    create_x_moment_mesh();
    // create_volume_constraint_mesh();
    create_flux_surface_meshes();
    create_flux_mesh();

    rebuild_global_mesh();

    cout << "setup_variable_and_function_pointers" << endl;
    setup_variable_and_function_pointers();

    cout << "set_boundary_conditions" << endl;
    set_boundary_conditions();

    // Setup equation numbering scheme
    cout << "Number of equations: " << endl;
    cout << assign_eqn_numbers() << endl;
    cout << "Number of unknowns: " << endl;
    cout << ndof() << endl;

    // DoubleVector residuals;
    // DenseDoubleMatrix jacobian;
    // DoubleVector residualsFD;
    // DenseDoubleMatrix jacobianFD(ndof());
    // DenseDoubleMatrix jacobian_error(ndof());

    // cout << "Get Jacobian" << endl;
    // get_jacobian(residuals, jacobian);
    // cout << "Get FD Jacobian" << endl;
    // get_fd_jacobian(residualsFD, jacobianFD);

    // for (unsigned i = 0; i < ndof(); i++)
    //{
    //  for (unsigned j = 0; j < ndof(); j++)
    //  {
    //    jacobian_error(i, j) = jacobian(i, j) - jacobianFD(i, j);

    //    if (abs(jacobian_error(i, j)) > 1e-3)
    //    {
    //      printf("i: %4u, j: %4u, act: %8.5f, exp: %8.5f \n",
    //             i,
    //             j,
    //             jacobian(i, j),
    //             jacobianFD(i, j));
    //    }
    //  }
    //}

    // for (unsigned i = 0; i < ndof(); i++)
    //{
    //  printf("i: %4u", i);
    //  for (unsigned j = ndof() - 4; j < ndof(); j++)
    //  {
    //    printf(", act: %9.2e", jacobian(j, i));
    //    // printf(", act: %9.2e, exp: %9.2e", jacobian(j, i), jacobianFD(j,
    //    // i));
    //  }
    //  cout << endl;
    //}

    // double j = 0;
    // cout << "Residuals" << endl;
    // for (unsigned i = 0; i < ndof(); i++)
    //{
    //  j = residuals[i];
    //  cout << "i: " << i << ", j: " << j << endl;
    //}
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::doc_fluid_mesh(string filename)
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

  // template<class ELEMENT>
  // void FingerProblem<ELEMENT>::doc_surface_mesh(string filename)
  //{
  //  ofstream output_stream;
  //  output_stream.open(filename);

  //  output_stream.close();
  //}

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::doc_volume(string filename)
  {
    ofstream output_stream;
    output_stream.open(filename, ofstream::app);

    output_stream << time() << ", ";
    output_stream << Volume_data_pt->value(0) << endl;

    output_stream.close();
  }

  // template<class ELEMENT>
  // void FingerProblem<ELEMENT>::doc_bubble_pressure(string filename)
  //{
  //  ofstream output_stream;
  //  output_stream.open(filename, ofstream::app);

  //  output_stream << time() << ", ";
  //  output_stream << Bubble_pressure_data_pt->value(0) << endl;

  //  output_stream.close();
  //}

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::doc_integrals(string filename)
  {
    ofstream output_stream;
    output_stream.open(filename, ofstream::app);

    output_stream << time() << ", ";
    output_stream << Outlet_integral_data_pt->value(0) << ", ";
    output_stream << Top_inlet_integral_data_pt->value(0) << ", ";
    output_stream << Bottom_inlet_integral_data_pt->value(0) << endl;

    output_stream.close();
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::compute_error_estimate(double& max_err,
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


  // template<class ELEMENT>
  // void FingerProblem<ELEMENT>::delete_surface_mesh()
  //{
  //  // How many surface elements are in the surface mesh
  //  unsigned n_element = Surface_mesh_pt->nelement();

  //  // Loop over the surface elements
  //  for (unsigned e = 0; e < n_element; e++)
  //  {
  //    // Delete surface element
  //    delete Surface_mesh_pt->element_pt(e);
  //  }

  //  // Wipe the mesh
  //  Surface_mesh_pt->flush_element_and_node_storage();
  //}

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::delete_volume_mesh()
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
  void FingerProblem<ELEMENT>::delete_x_moment_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = X_moment_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete X_moment_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    X_moment_mesh_pt->flush_element_and_node_storage();
  }

  // template<class ELEMENT>
  // void FingerProblem<ELEMENT>::delete_volume_constraint_mesh()
  //{
  //  // How many elements are in the volume mesh
  //  unsigned n_element = Volume_constraint_mesh_pt->nelement();

  //  // Loop over the elements
  //  for (unsigned e = 0; e < n_element; e++)
  //  {
  //    // Delete element
  //    delete Volume_constraint_mesh_pt->element_pt(e);
  //  }

  //  // Wipe the mesh
  //  Volume_constraint_mesh_pt->flush_element_and_node_storage();
  //}

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::delete_flux_surface_meshes()
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

    // How many elements are in the volume mesh
    n_element = Top_inlet_surface_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete Top_inlet_surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Top_inlet_surface_mesh_pt->flush_element_and_node_storage();

    // How many elements are in the volume mesh
    n_element = Bottom_inlet_surface_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Delete element
      delete Bottom_inlet_surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Bottom_inlet_surface_mesh_pt->flush_element_and_node_storage();
  }

  template<class ELEMENT>
  void FingerProblem<ELEMENT>::delete_flux_mesh()
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
  void FingerProblem<ELEMENT>::delete_mesh_pt(Mesh* mesh_pt)
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
