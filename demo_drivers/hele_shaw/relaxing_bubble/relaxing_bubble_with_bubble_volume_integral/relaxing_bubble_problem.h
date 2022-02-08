#ifndef OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER
#define OOMPH_RELAXING_BUBBLE_PROBLEM_HEADER

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "solid.h"
#include "relaxing_bubble_parameters.h"
#include "compare_nodes.h"

template<class ELEMENT, class BUBBLE_ELEMENT>
class RelaxingBubbleProblem : public Problem
{
  /// Variables and enums
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
  /// Polygons
  TriangleMeshPolygon* Outer_boundary_polygon_pt;
  TriangleMeshPolygon* Surface_polygon_pt;

  /// Meshes
  RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;
  RefineableSolidTriangleMesh<BUBBLE_ELEMENT>* Bubble_mesh_pt;
  Mesh* Surface_mesh_pt;
  Mesh* Info_mesh_pt;
  Mesh* Volume_constraint_mesh_pt;

  /// Data
  Data* Integral_data_pt;
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

  /// Compare given Jacobian to a numerical one
  void debug_jacobian()
  {
    DoubleVector residuals;
    DenseDoubleMatrix jacobian;
    DoubleVector residualsFD;
    DenseDoubleMatrix jacobianFD(ndof());

    cout << "Get Jacobian" << endl;
    get_jacobian(residuals, jacobian);
    get_fd_jacobian(residualsFD, jacobianFD);

    double err;
    for (unsigned i = 0; i < ndof(); i++)
    {
      for (unsigned j = 0; j < ndof(); j++)
      {
        err = jacobian(j, i) - jacobianFD(j, i);
        if (abs(err) > 1e-4)
        {
          printf("i: %4u, j: %4u", i, j);
          printf(
            ", act: %10.5g, exp: %10.5g", jacobian(j, i), jacobianFD(j, i));
          printf(", dif: %10.5g", err);
          cout << endl;
        }
      }
    }

    cout << "Final rows" << endl;
    for (unsigned i = 0; i < ndof(); i++)
    {
      for (unsigned j = ndof() - 5; j < ndof(); j++)
      {
        printf("%10.5g, ", jacobian(j, i));
      }
      printf("\n");
    }

    cout << "Final columns" << endl;
    for (unsigned j = 0; j < ndof(); j++)
    {
      for (unsigned i = ndof() - 5; i < ndof(); i++)
      {
        printf("%10.5g, ", jacobian(j, i));
      }
      printf("\n");
    }
  }

  void print_cross_mesh_jacobian_terms()
  {
    DoubleVector residuals;
    DenseDoubleMatrix jacobian;

    cout << "Get Jacobian" << endl;
    get_jacobian(residuals, jacobian);

    cout << "Bubble rows" << endl;
    for (unsigned i = 0; i < 1277; i++)
    {
      for (unsigned j = 1277; j < 1277+366; j++)
      {
        printf("%1u", abs(jacobian(j, i))>0.0);
      }
      printf("\n");
    }

  }

private:
  void create_data_pointers();

  void generate_mesh();
  void generate_outer_boundary_polygon();
  void generate_surface_polygon();
  void generate_fluid_mesh();
  void generate_surface_mesh();

  void generate_bubble_mesh();
  void get_bubble_bounding_nodes(Vector<Vector<Node*>>& boundary_nodes_pt);
  void sort_bubble_bounding_nodes(Vector<Vector<Node*>>& boundary_nodes_pt);
  void point_bubble_boundary_nodes_at_fluid_ones(
    Vector<Vector<Node*>>& fluid_boundary_nodes_pt);

  void set_variable_and_function_pointers();

  void set_boundary_conditions();
  void fill_in_bubble_boundary_map(map<unsigned, bool>& is_on_bubble_bound);

  void actions_before_adapt();
  void actions_after_adapt();

  void doc_fluid_mesh(string filename);
  void doc_bubble_mesh(string filename);
  void doc_surface_mesh(string filename);
  void doc_integral(string filename);

  void compute_error_estimate(double& max_err, double& min_err);

  void delete_surface_mesh();
};

template<class ELEMENT, class BUBBLE_ELEMENT>
RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::RelaxingBubbleProblem()
{
  bool adaptive_timestepping = false;
  add_time_stepper_pt(new BDF<2>(adaptive_timestepping));

  create_data_pointers();

  generate_mesh();

  set_variable_and_function_pointers();

  set_boundary_conditions();

  // Setup equation numbering scheme
  cout << "Number of equations: " << endl;
  cout << assign_eqn_numbers() << endl;
  cout << "Number of unknowns: " << endl;
  cout << ndof() << endl;
}

/// Use a newton solve to set the initial conditions
template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::
  solve_for_initial_conditions(DocInfo& doc_info)
{
  // cout<<"Steady Newton solve"<<endl;
  // const unsigned max_adapt = 0;
  // steady_newton_solve(max_adapt);
  // newton_solve();
  double dt = 5e-3;
  initialise_dt(dt);
  assign_initial_values_impulsive(dt);
  Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

  doc_solution(doc_info);
};

/// Iterate forward in time
template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::iterate_timestepper(
  const double& t_step, const double& t_final, DocInfo& doc_info)
{
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
    unsteady_newton_solve(t_step, max_adapt, is_first_step);
    if (is_first_step)
    {
      is_first_step = false;
    }
    Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

    doc_solution(doc_info);
  }
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::create_data_pointers()
{
  cout << "Create data" << endl;
  Integral_data_pt = new Data(1);
  Bubble_pressure_data_pt = new Data(1);
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT,
                           BUBBLE_ELEMENT>::generate_outer_boundary_polygon()
{
  cout << "Generate outer boundary polygon" << endl;
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

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::generate_surface_polygon()
{
  cout << "Generate surface polygon" << endl;
  double x_center = 0.5;
  double y_center = 0.5;
  double major_radius = 0.4;
  double minor_radius = pow(0.3, 2.0) / major_radius;
  unsigned npoints = 16;
  double zeta_step = MathematicalConstants::Pi / double(npoints - 1);

  // Intrinsic coordinate along GeomObject defining the bubble
  Vector<double> zeta(1);

  // Position vector to GeomObject defining the bubble
  Vector<double> coord(2);

  Ellipse* ellipse_pt = new Ellipse(major_radius, minor_radius);

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

  surface_polyline_pt[0]->set_maximum_length(0.02);
  surface_polyline_pt[1]->set_maximum_length(0.02);

  // Create closed polygon from two polylines
  Surface_polygon_pt =
    new TriangleMeshPolygon(surface_polyline_pt, bubble_center);

  Surface_polygon_pt->set_polyline_refinement_tolerance(0.05);
  Surface_polygon_pt->set_polyline_unrefinement_tolerance(0.01);
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::generate_mesh()
{
  generate_outer_boundary_polygon();
  generate_surface_polygon();

  cout << "Generate fluid mesh" << endl;
  generate_fluid_mesh();
  cout << "Generate bubble mesh" << endl;
  generate_bubble_mesh();

  //  DocInfo doc_info;
  //  doc_info.set_directory("TEST/");
  //  doc_solution(doc_info);
  //
  cout << "Generate surface mesh" << endl;

  this->Info_mesh_pt = new Mesh;
  this->Info_mesh_pt->add_element_pt(new InfoElement(Integral_data_pt));

  Volume_constraint_mesh_pt = new Mesh;
  MyConstraintElement* vol_constraint_element =
    new MyConstraintElement(relaxing_bubble::target_bubble_volume_pt,
                            Integral_data_pt->value_pt(0),
                            Bubble_pressure_data_pt,
                            0);
  Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);

  Surface_mesh_pt = new Mesh;
  generate_surface_mesh();

  add_sub_mesh(Fluid_mesh_pt);
  add_sub_mesh(Bubble_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);
  add_sub_mesh(Info_mesh_pt);
  add_sub_mesh(Volume_constraint_mesh_pt);

  build_global_mesh();

  Fluid_mesh_pt->output("fluid.dat");
  Fluid_mesh_pt->output_boundaries("fluid-boundaries.dat");
  cout << "Number of fluid mesh elements: " << Fluid_mesh_pt->nelement()
       << endl;
  Bubble_mesh_pt->output("bubble.dat");
  Bubble_mesh_pt->output_boundaries("bubble-boundaries.dat");
  cout << "Number of bubble mesh elements: " << Bubble_mesh_pt->nelement()
       << endl;
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::generate_fluid_mesh()
{
  // Target area for initial mesh
  double uniform_element_area = 5e-2;

  TriangleMeshClosedCurve* outer_closed_curve_pt = Outer_boundary_polygon_pt;
  Vector<TriangleMeshClosedCurve*> surface_closed_curve_pt(1);
  surface_closed_curve_pt[0] = Surface_polygon_pt;

  // Use the TriangleMeshParameters object for gathering all
  // the necessary arguments for the TriangleMesh object
  TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

  // Define the holes on the boundary
  triangle_mesh_parameters.internal_closed_curve_pt() = surface_closed_curve_pt;

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

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::generate_surface_mesh()
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

  unsigned i_element = 0;
  MyConstraintElement* vol_constraint_element =
    dynamic_cast<MyConstraintElement*>(
      Volume_constraint_mesh_pt->element_pt(i_element));

  bool fd_jacobian = true;
  vol_constraint_element->add_external_data(Integral_data_pt, fd_jacobian);
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::generate_bubble_mesh()
{
  Vector<Vector<Node*>> boundary_nodes_pt(2);

  cout << "Get bubble bounding nodes" << endl;
  get_bubble_bounding_nodes(boundary_nodes_pt);

  cout << "Sort bubble bounding nodes" << endl;
  sort_bubble_bounding_nodes(boundary_nodes_pt);

  cout << "Create mesh" << endl;
  /// Triangle requires at least two polylines for a single polygon
  Vector<TriangleMeshCurveSection*> bubble_polyline_pt(2);

  // Vertex coordinates
  Vector<Vector<double>> bubble_vertex;

  const unsigned boundary_id[2] = {First_bubble_boundary_id,
                                   Second_bubble_boundary_id};

  for (unsigned i = 0; i < 2; i++)
  {
    unsigned npoints = boundary_nodes_pt[i].size();

    bubble_vertex.resize(npoints);

    for (unsigned j = 0; j < npoints; j++)
    {
      bubble_vertex[j].resize(2);
      boundary_nodes_pt[i][j]->position(bubble_vertex[j]);
      cout << "(x,y) = (" << bubble_vertex[j][0] << ", " << bubble_vertex[j][1]
           << ") " << endl;
    }

    // Build the 1st bubble polyline
    bubble_polyline_pt[i] =
      new TriangleMeshPolyLine(bubble_vertex, boundary_id[i]);
  }

  // Define coordinates of a point inside the bubble
  // TODO Double check
  double x_center = 0.5;
  double y_center = 0.5;
  Vector<double> bubble_center(2);
  bubble_center[0] = x_center;
  bubble_center[1] = y_center;

  // Create closed polygon from two polylines
  TriangleMeshPolygon* bubble_polygon_pt =
    new TriangleMeshPolygon(bubble_polyline_pt, bubble_center);

  TriangleMeshClosedCurve* bubble_closed_curve_pt = bubble_polygon_pt;

  // Use the TriangleMeshParameters object for gathering all
  // the necessary arguments for the TriangleMesh object
  TriangleMeshParameters triangle_mesh_parameters(bubble_closed_curve_pt);

  // Disable creation of vertices between the given boundary nodes
  triangle_mesh_parameters
    .disable_automatic_creation_of_vertices_on_boundaries();

  // Target area for initial mesh
  double uniform_element_area = 5e-2;

  // Define the maximum element areas
  triangle_mesh_parameters.element_area() = uniform_element_area;

  cout << "Call RefineableSolidTriangleMesh" << endl;
  // Create the mesh
  Bubble_mesh_pt = new RefineableSolidTriangleMesh<BUBBLE_ELEMENT>(
    triangle_mesh_parameters, this->time_stepper_pt());
  Bubble_mesh_pt->disable_projection();

  // Set error estimator for fluid mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Bubble_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  Bubble_mesh_pt->max_permitted_error() = 1e-3;
  Bubble_mesh_pt->min_permitted_error() = 5e-6;
  Bubble_mesh_pt->max_element_size() = 5e-2;
  Bubble_mesh_pt->min_element_size() = 1e-4;

  cout << "Point bubble boundary nodes at fluid ones" << endl;
  point_bubble_boundary_nodes_at_fluid_ones(boundary_nodes_pt);
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::get_bubble_bounding_nodes(
  Vector<Vector<Node*>>& boundary_nodes_pt)
{
  const unsigned boundary_id[2] = {First_bubble_boundary_id,
                                   Second_bubble_boundary_id};

  // Collect vectors of the nodes on the mesh boundaries
  for (unsigned i = 0; i < 2; i++)
  {
    unsigned b = boundary_id[i];

    unsigned n_node = Fluid_mesh_pt->nboundary_node(b);

    boundary_nodes_pt[i].resize(n_node);
    for (unsigned n = 0; n < n_node; n++)
    {
      boundary_nodes_pt[i][n] = Fluid_mesh_pt->boundary_node_pt(b, n);
    }
  }
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::sort_bubble_bounding_nodes(
  Vector<Vector<Node*>>& boundary_nodes_pt)
{
  std::sort(boundary_nodes_pt[0].begin(),
            boundary_nodes_pt[0].end(),
            CompareNodeCoordinatesX());

  std::sort(boundary_nodes_pt[1].rbegin(),
            boundary_nodes_pt[1].rend(),
            CompareNodeCoordinatesX());

  // Get rid of extra nodes!
  for (unsigned i = 0; i < 2; i++)
  {
    unsigned npoints = boundary_nodes_pt[i].size();
    for (unsigned j = 0; j < floor(0.5 * (npoints - 1)); j++)
    {
      boundary_nodes_pt[i].erase(boundary_nodes_pt[i].begin() + 1 + j);
    }
  }
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::
  point_bubble_boundary_nodes_at_fluid_ones(
    Vector<Vector<Node*>>& fluid_boundary_nodes_pt)
{
  const unsigned boundary_id[2] = {First_bubble_boundary_id,
                                   Second_bubble_boundary_id};

  // Collect vectors of the nodes on the mesh boundaries
  for (unsigned i = 0; i < 2; i++)
  {
    unsigned b = boundary_id[i];

    unsigned n_node = Bubble_mesh_pt->nboundary_node(b);
    for (unsigned n = 0; n < n_node; n++)
    {
      dynamic_cast<Mesh*>(Bubble_mesh_pt)->boundary_node_pt(b, n) =
        fluid_boundary_nodes_pt[i][n];
    }
  }
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT,
                           BUBBLE_ELEMENT>::set_variable_and_function_pointers()
{
  Integral_data_pt->set_value(0, (*relaxing_bubble::target_bubble_volume_pt));

  Bubble_pressure_data_pt->set_value(0, 16.0 / 3.0);
  relaxing_bubble::bubble_pressure_pt = Bubble_pressure_data_pt->value_pt(0);


  bool fd_jacobian = true;
  /// Set fluid mesh function pointers
  unsigned n_element = Fluid_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

    // Set the constitutive law for pseudo-elastic mesh deformation
    el_pt->constitutive_law_pt() = relaxing_bubble::constitutive_law_pt;
    el_pt->upper_wall_fct_pt() = relaxing_bubble::upper_wall_fct;
    el_pt->add_external_data(Integral_data_pt, fd_jacobian);
    el_pt->add_external_data(Bubble_pressure_data_pt, fd_jacobian);
  }

  /// Set fluid mesh function pointers
  n_element = Bubble_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    BUBBLE_ELEMENT* el_pt =
      dynamic_cast<BUBBLE_ELEMENT*>(Bubble_mesh_pt->element_pt(e));

    // Set the constitutive law for pseudo-elastic mesh deformation
    el_pt->constitutive_law_pt() = relaxing_bubble::constitutive_law_pt;
    el_pt->integrand_fct_pt() = relaxing_bubble::depth_fct;
    el_pt->set_external_data_pt(Integral_data_pt);
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

    interface_element_pt->add_external_data(Integral_data_pt, fd_jacobian);
    interface_element_pt->add_external_data(Bubble_pressure_data_pt,
                                            fd_jacobian);
  }

  MyConstraintElement* vol_constraint_element =
    dynamic_cast<MyConstraintElement*>(
      Volume_constraint_mesh_pt->element_pt(0));
  vol_constraint_element->add_external_data(Integral_data_pt, fd_jacobian);
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::set_boundary_conditions()
{
  cout << "Set boundary conditions" << endl;

  Integral_data_pt->unpin(0);
  Bubble_pressure_data_pt->unpin(0);

  Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
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
        solid_node_pt->pin_position(0);
        solid_node_pt->pin_position(1);
      }
    }
  } // end loop over boundaries

  // const unsigned boundary_id[2] = {First_bubble_boundary_id,
  //                                 Second_bubble_boundary_id};
  // nbound = Bubble_mesh_pt->nboundary();
  // for (unsigned ibound = 0; ibound < nbound; ibound++)
  //{
  //  unsigned num_nod = Bubble_mesh_pt->nboundary_node(boundary_id[ibound]);
  //  for (unsigned inod = 0; inod < num_nod; inod++)
  //  {
  //    // Get node
  //    Node* node_pt = Bubble_mesh_pt->boundary_node_pt(boundary_id[ibound],
  //    inod);

  //    // Pin pseudo-solid positions apart from bubble boundary which
  //    // we allow to move
  //    SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(node_pt);

  //    solid_node_pt->unpin_position(0);
  //    solid_node_pt->unpin_position(1);
  //  }
  //}

  /// Single point Dirichlet boundary condition
  Node* node_pt = Fluid_mesh_pt->boundary_node_pt(0, 0);
  const double fixed_pressure = 0.0;
  node_pt->set_value(0, fixed_pressure);
  node_pt->pin(0);


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

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::
  fill_in_bubble_boundary_map(map<unsigned, bool>& is_on_bubble_bound)
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

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::doc_solution(
  DocInfo& doc_info)
{
  string doc_directory = doc_info.directory();

  doc_fluid_mesh(doc_directory + "soln" + to_string(doc_info.number()) +
                 ".dat");
  // doc_surface_mesh(doc_directory + "surface" + doc_info.number() + ".dat");

  doc_bubble_mesh(doc_directory + "bubble" + to_string(doc_info.number()) +
                  ".dat");

  doc_integral(doc_directory + "integral.dat");

  doc_info.number()++;
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::actions_before_adapt()
{
  delete_surface_mesh();

  delete this->Info_mesh_pt->element_pt(0);
  Info_mesh_pt->flush_element_and_node_storage();

  rebuild_global_mesh();
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::actions_after_adapt()
{
  Integral_data_pt = new Data(1);

  this->Info_mesh_pt->add_element_pt(new InfoElement(Integral_data_pt));

  generate_surface_mesh();

  rebuild_global_mesh();

  set_variable_and_function_pointers();

  set_boundary_conditions();

  // Setup equation numbering scheme
  cout << "Number of equations: " << endl;
  cout << assign_eqn_numbers() << endl;
  cout << "Number of unknowns: " << endl;
  cout << ndof() << endl;
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::doc_fluid_mesh(
  string filename)
{
  double max_err;
  double min_err;
  compute_error_estimate(max_err, min_err);

  ofstream output_stream;
  output_stream.open(filename);
  unsigned npoints = 1;
  Fluid_mesh_pt->output(output_stream, npoints);
  output_stream.close();
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::doc_bubble_mesh(
  string filename)
{
  // double max_err;
  // double min_err;
  // compute_error_estimate(max_err, min_err);

  ofstream output_stream;
  output_stream.open(filename);
  unsigned npoints = 1;
  Bubble_mesh_pt->output(output_stream, npoints);
  output_stream.close();
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::doc_integral(
  string filename)
{
  ofstream output_stream;
  output_stream.open(filename, ios_base::app);
  output_stream << "t: " << time() << ", I: " << Integral_data_pt->value(0)
                << endl;
  output_stream.close();
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::doc_surface_mesh(
  string filename)
{
  ofstream output_stream;
  output_stream.open(filename);

  output_stream.close();
}

template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::compute_error_estimate(
  double& max_err, double& min_err)
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


template<class ELEMENT, class BUBBLE_ELEMENT>
void RelaxingBubbleProblem<ELEMENT, BUBBLE_ELEMENT>::delete_surface_mesh()
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

#endif
