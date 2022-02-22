#ifndef OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER
#define OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "info_element.h"

#include "inject_air_parameters.h"
#include "my_constraint_elements.h"
#include "modified_volume_constraint_elements_with_integrals.h"

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
    Mesh* Volume_constraint_mesh_pt;

    VolumeConstraintElement* Vol_constraint_el_pt;

    Data* Bubble_pressure_data_pt;

    /// Functions
  public:
    /// Constructor
    RelaxingBubbleProblem();

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

    void actions_before_implicit_timestep()
    {
      inject_air::volume_fct(time(), *inject_air::target_bubble_volume_pt);
    }

    void actions_before_newton_solve()
    {
      // debug_jacobian();
    }

    void actions_after_newton_solve()
    {
      // debug_jacobian();
    }


    void debug_jacobian()
    {
      DoubleVector residuals;
      DenseDoubleMatrix jacobian;
      DoubleVector residualsFD;
      DenseDoubleMatrix jacobianFD(ndof());

      cout << "Get Jacobian" << endl;
      get_jacobian(residuals, jacobian);
      get_fd_jacobian(residualsFD, jacobianFD);

      for (unsigned i = 0; i < ndof(); i++)
      {
        for (unsigned j = ndof() - 3; j < ndof(); j++)
        {
          if (abs(jacobian(j, i) - jacobianFD(j, i)) > 1e-6)
          {
            printf("i: %4u, j: %4u, act: %8.5f, exp: %8.5f\n",
                   i,
                   j,
                   jacobian(j, i),
                   jacobianFD(j, i));
          }
        }
      }

      cout << "Residuals" << endl;
      for (unsigned i = 0; i < ndof(); i++)
      {
        cout << "i: " << i << ", res: " << residuals[i] << endl;
      }
    }

  private:
    void create_mesh();
    void create_data();
    void create_outer_boundary_polygon();
    void create_surface_polygon();
    void create_fluid_mesh();
    void create_surface_mesh();
    void create_volume_constraint_mesh();

    void set_variable_and_function_pointers();

    void set_boundary_conditions();
    void fill_in_bubble_boundary_map(map<unsigned, bool>& is_on_bubble_bound);

    void actions_before_adapt();
    void actions_after_adapt();

    void doc_fluid_mesh(string filename);
    void doc_surface_mesh(string filename);
    void doc_bubble_pressure(string filename);
    void doc_boundary(string filename);

    void compute_error_estimate(double& max_err, double& min_err);

    void delete_surface_mesh();
    void delete_volume_constraint_mesh();
    void delete_mesh_pt(Mesh* mesh_pt);
  };

  template<class ELEMENT>
  RelaxingBubbleProblem<ELEMENT>::RelaxingBubbleProblem()
  {
    cout << "Constructor" << endl;
    bool adaptive_timestepping = false;
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

    double dt = 1e-2;
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
    unsigned adapt_interval;
    if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
    {
      adapt_interval = 4;
    }
    else
    {
      adapt_interval = 1;
    }
    bool remesh_initial_condition = false;
    for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++)
    {
      cout << "t: " << time() << endl;

      // max_newton_iterations() = 1;
      // newton_solver_tolerance() = 1e10;
      if (i_timestep % adapt_interval == 0 && !is_first_step)
      {
        max_adapt = 1;
      }
      else
      {
        max_adapt = 0;
      }


      /// Allow (significant??) adapting on the first step

      if (is_first_step && remesh_initial_condition)
      {
        max_adapt = 1;
      }

      // this->linear_solver_pt() = new FD_LU;

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
    const double x_center = 0.5;
    const double y_center = 0.5;
    double minor_radius =
      (-*inject_air::target_bubble_volume_pt) /
      (MathematicalConstants::Pi * inject_air::major_radius);
    unsigned npoints = 32;
    double zeta_step = MathematicalConstants::Pi / double(npoints - 1);

    // Intrinsic coordinate along GeomObject defining the bubble
    Vector<double> zeta(1);

    // Position vector to GeomObject defining the bubble
    Vector<double> coord(2);

    Ellipse* ellipse_pt = new Ellipse(inject_air::major_radius, minor_radius);

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

    double max_length = 0.02;
    surface_polyline_pt[0]->set_maximum_length(max_length);
    surface_polyline_pt[1]->set_maximum_length(max_length);

    // Create closed polygon from two polylines
    Surface_polygon_pt =
      new TriangleMeshPolygon(surface_polyline_pt, bubble_center);

    Surface_polygon_pt->set_polyline_refinement_tolerance(0.08);
    Surface_polygon_pt->set_polyline_unrefinement_tolerance(0.04);
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_data()
  {
    Bubble_pressure_data_pt = new Data(1);
    Bubble_pressure_data_pt->set_value(0, 0.0);
    add_global_data(Bubble_pressure_data_pt);
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_mesh()
  {
    create_fluid_mesh();

    Surface_mesh_pt = new Mesh;
    create_surface_mesh();

    Volume_constraint_mesh_pt = new Mesh;

    unsigned index_of_traded_pressure = 0;
    Vol_constraint_el_pt =
      new VolumeConstraintElement(inject_air::target_bubble_volume_pt,
                                  Bubble_pressure_data_pt,
                                  index_of_traded_pressure);

    // Add volume constraint element to the mesh
    Volume_constraint_mesh_pt->add_element_pt(Vol_constraint_el_pt);

    create_volume_constraint_mesh();

    add_sub_mesh(Fluid_mesh_pt);
    add_sub_mesh(Surface_mesh_pt);
    add_sub_mesh(Volume_constraint_mesh_pt);

    build_global_mesh();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::create_fluid_mesh()
  {
    // Target area for initial mesh
    double uniform_element_area = 1e-2;


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
    if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
    {
      Fluid_mesh_pt->max_permitted_error() = 1e-2;
      Fluid_mesh_pt->min_element_size() = 1e-2;
    }
    else
    {
      Fluid_mesh_pt->max_permitted_error() = 1e-2;
      // Fluid_mesh_pt->min_permitted_error() = 1e-6;
      // Fluid_mesh_pt->max_element_size() = 5e-1;
      Fluid_mesh_pt->min_element_size() = 1e-3;
    }
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
  void RelaxingBubbleProblem<ELEMENT>::create_volume_constraint_mesh()
  {
    // Loop over the free surface boundaries
    // unsigned nb=Fluid_mesh_pt->nboundary();
    for (unsigned b = First_bubble_boundary_id;
         b < Second_bubble_boundary_id + 1;
         b++)
    {
      std::cout << "Setting constraint elements on boundary " << b << std::endl;
      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element = Fluid_mesh_pt->nboundary_element(b);

      // Loop over the bulk fluid elements adjacent to boundary b?
      for (unsigned e = 0; e < n_element; e++)
      {
        // Get pointer to the bulk fluid element that is
        // adjacent to boundary b
        ELEMENT* bulk_elem_pt =
          dynamic_cast<ELEMENT*>(Fluid_mesh_pt->boundary_element_pt(b, e));

        // Find the index of the face of element e along boundary b
        int face_index = Fluid_mesh_pt->face_index_at_boundary(b, e);

        // Create new element
        HeleShawVolumeConstraintElement<ELEMENT>* el_pt =
          new HeleShawVolumeConstraintElement<ELEMENT>(bulk_elem_pt,
                                                       face_index);

        // Set the "master" volume constraint element
        el_pt->set_volume_constraint_element(Vol_constraint_el_pt);

        // Add it to the mesh
        Volume_constraint_mesh_pt->add_element_pt(el_pt);
      }
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::set_variable_and_function_pointers()
  {
    bool fd_jacobian = true;

    inject_air::bubble_pressure_pt = Bubble_pressure_data_pt->value_pt(0);

    /// Set fluid mesh function pointers
    unsigned n_element = Fluid_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast from GeneralisedElement to the present element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

      // Set the constitutive law for pseudo-elastic mesh deformation
      el_pt->constitutive_law_pt() = inject_air::constitutive_law_pt;
      el_pt->upper_wall_fct_pt() = inject_air::upper_wall_fct;
      el_pt->add_external_data(Bubble_pressure_data_pt, fd_jacobian);
    }

    n_element = Surface_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      HeleShawInterfaceElement<ELEMENT>* interface_element_pt =
        dynamic_cast<HeleShawInterfaceElement<ELEMENT>*>(
          Surface_mesh_pt->element_pt(e));
      interface_element_pt->ca_inv_pt() = &inject_air::ca_inv;
      interface_element_pt->st_pt() = &inject_air::st;
      interface_element_pt->aspect_ratio_pt() = &inject_air::alpha;
      interface_element_pt->upper_wall_fct_pt() = inject_air::upper_wall_fct;
      interface_element_pt->wall_speed_fct_pt() = inject_air::wall_speed_fct;
      interface_element_pt->bubble_pressure_fct_pt() =
        inject_air::bubble_pressure_fct;

      interface_element_pt->add_external_data(Bubble_pressure_data_pt,
                                              fd_jacobian);
    }
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::set_boundary_conditions()
  {
    Bubble_pressure_data_pt->unpin(0);

    // Map to record if a given boundary is on a bubble or not
    map<unsigned, bool> is_on_bubble_bound;
    fill_in_bubble_boundary_map(is_on_bubble_bound);

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
          node_pt->set_value(0, 0.0);
          node_pt->pin(0);

          solid_node_pt->pin_position(0);
          solid_node_pt->pin_position(1);
        }
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

    doc_bubble_pressure(doc_directory + "bubble_pressure.dat");

    doc_info.number()++;
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::actions_before_adapt()
  {
    cout << "Actions before adapt" << endl;

    delete_surface_mesh();
    delete_volume_constraint_mesh();

    flush_global_data();

    rebuild_global_mesh();
  }

  template<class ELEMENT>
  void RelaxingBubbleProblem<ELEMENT>::actions_after_adapt()
  {
    cout << "Actions after adapt" << endl;

    create_surface_mesh();
    create_volume_constraint_mesh();

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
  void RelaxingBubbleProblem<ELEMENT>::doc_bubble_pressure(string filename)
  {
    ofstream output_stream;
    output_stream.open(filename, ofstream::app);

    output_stream << time() << ", ";
    output_stream << Bubble_pressure_data_pt->value(0) << endl;

    output_stream.close();
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
  void RelaxingBubbleProblem<ELEMENT>::delete_volume_constraint_mesh()
  {
    // How many elements are in the volume mesh
    unsigned n_element = Volume_constraint_mesh_pt->nelement();

    // Loop over the elements
    for (unsigned e = 1; e < n_element; e++)
    {
      // Delete element
      delete Volume_constraint_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Volume_constraint_mesh_pt->flush_element_and_node_storage();
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
