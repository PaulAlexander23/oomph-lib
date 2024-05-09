//==start_of_problem_class============================================
/// Problem class to simulate inviscid bubble propagating along 2D channel
//====================================================================
template<class ELEMENT>
class BubbleInChannelProblem : public Problem
{
  /// Make a copy for using in bifurcation tracking


public:
  double Q_jack;

  /// Constructor
  BubbleInChannelProblem();


  Problem* make_copy()
  {
    // Make a copy based on the current parameters
    return (new BubbleInChannelProblem());
  }


  /// Destructor
  ~BubbleInChannelProblem()
  {
    // Fluid timestepper
    delete this->time_stepper_pt(0);

    // Kill data associated with outer boundary
    unsigned n = Outer_boundary_polyline_pt->npolyline();
    for (unsigned j = 0; j < n; j++)
    {
      delete Outer_boundary_polyline_pt->polyline_pt(j);
    }
    delete Outer_boundary_polyline_pt;

    // Kill data associated with bubbles
    unsigned n_bubble = Bubble_polygon_pt.size();
    for (unsigned ibubble = 0; ibubble < n_bubble; ibubble++)
    {
      unsigned n = Bubble_polygon_pt[ibubble]->npolyline();
      for (unsigned j = 0; j < n; j++)
      {
        delete Bubble_polygon_pt[ibubble]->polyline_pt(j);
      }
      delete Bubble_polygon_pt[ibubble];
    }


    // Flush element of free surface elements
    delete_free_surface_elements();
    delete Free_surface_mesh_pt;
    delete_volume_constraint_elements();
    delete Volume_constraint_mesh_pt;

    delete_CoM_X_constraint_elements();
    delete_CoM_Y_constraint_elements();
    delete_inflow_elements();

    delete CoM_X_constraint_mesh_pt;
    delete CoM_Y_constraint_mesh_pt;
    delete Inflow_mesh_pt;

    // Delete error estimator
    delete Fluid_mesh_pt->spatial_error_estimator_pt();

    // Delete fluid mesh
    delete Fluid_mesh_pt;

    // Delete the global pressure bubble data
    delete Bubble_pressure_data_pt;

    // Kill const eqn
    delete Problem_Parameter::Constitutive_law_pt;
  }


  /// Actions before adapt: Wipe the mesh of free surface elements
  void actions_before_adapt()
  {
    // Kill the  elements and wipe surface mesh
    delete_free_surface_elements();
    delete_inflow_elements();
    delete_volume_constraint_elements();
    delete_CoM_X_constraint_elements();
    delete_CoM_Y_constraint_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();

  } // end of actions_before_adapt


  /// Actions after adapt: Rebuild the mesh of free surface elements
  void actions_after_adapt()
  {
    // Create the elements that impose the displacement constraint
    create_free_surface_elements();
    create_inflow_elements();
    create_volume_constraint_elements();
    create_CoM_X_constraint_elements();
    create_CoM_Y_constraint_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();


    complete_problem_setup();

    pin_tangential_lagrange();
    redo_equation_numbering();
  } // end of actions_after_adapt


  double global_temporal_error_norm()
  {
    std::cout << "Calling global error norm" << std::endl;
    double global_error = 0.0;

    unsigned count = 0;
    unsigned nbound = Fluid_mesh_pt->nboundary();
    for (unsigned ibound = 0; ibound < nbound; ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);

        // Get temporal error in position from x and y coordinate of
        // each node on the boundary.
        if (1)
        {
          double errx =
            nod_pt->position_time_stepper_pt()->temporal_error_in_position(
              nod_pt, 0);
          double erry =
            nod_pt->position_time_stepper_pt()->temporal_error_in_position(
              nod_pt, 1);

          // Add the square of the individual error to the global error
          count++;
          global_error += errx * errx + erry * erry;
        }
      }
    }
    oomph_info << "Global temporal error norm: " << global_error << " over "
               << count << " interface nodes" << std::endl;
    // Divide by the number of nodes
    global_error /= double(count);

    // Return square root...
    return std::sqrt(global_error);
  }
  /// Update the after solve (empty)
  void actions_after_newton_solve() {}

  //    void actions_before_newton_convergence_check()
  //{
  //      std::cout << "Q is " << get_Q() << "<n";
  //}

  /// Update the problem specs before solve
  void actions_before_newton_solve() {}

  /// \short Set boundary conditions and complete the build of all elements
  void complete_problem_setup();

  // Custom function to perform continuation in 1 parameter
  void jack_solve(double Q, bool jack_arc, unsigned n_iter, double sign);

  // Custom function to perform continuation in 2 parameters
  void jack_solve(double Q, double V, bool jack_arc, unsigned n_iter, int sign);

  void dump_it(ofstream& dump_file);

  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  /// Compute the error estimates and assign to elements for plotting
  void compute_error_estimate(double& max_err, double& min_err);

  void reset_lagrangian_coordinates()
  {
    Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  }

  double get_bubble_pressure()
  {
    return Bubble_pressure_data_pt->value(0);
  }

  void create_inflow_elements()
  {
    unsigned boundary_for_flux_elements = Inflow_boundary_id;
    unsigned n_inflow_element =
      Fluid_mesh_pt->nboundary_element(boundary_for_flux_elements);
    for (unsigned e = 0; e < n_inflow_element; e++)
    {
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(boundary_for_flux_elements, e));
      // Find the index of the face of element e along boundary b
      int face_index =
        Fluid_mesh_pt->face_index_at_boundary(boundary_for_flux_elements, e);
      HeleShawFluxElement<ELEMENT>* el_pt =
        new HeleShawFluxElement<ELEMENT>(bulk_elem_pt, face_index);
      Inflow_mesh_pt->add_element_pt(el_pt);
      el_pt->flux_fct_pt() = &Problem_Parameter::normal_flux_behind_bubble;

      /// This one is important!
      el_pt->add_external_data(G_data_pt, true);
    }

    cout << "Poisson Flux Elements created on main inflow boundary" << endl;
  }


  void delete_inflow_elements()
  {
    unsigned n_element = Inflow_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      delete Inflow_mesh_pt->element_pt(e);
    }
    Inflow_mesh_pt->flush_element_and_node_storage();
  }
  void delete_outflow_elements()
  {
    unsigned n_element = Outflow_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      delete Outflow_mesh_pt->element_pt(e);
    }
    Outflow_mesh_pt->flush_element_and_node_storage();
  }


  void pin_all_nodal_positions()
  {
    unsigned n_node = Fluid_mesh_pt->nnode();
    for (unsigned inod = 0; inod < n_node; inod++)
    {
      Node* nod_pt = Fluid_mesh_pt->node_pt(inod);
      SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);

      solid_node_pt->pin_position(0);
      solid_node_pt->pin_position(1);
    }
    std::cout << "Pin nodal positions " << std::endl;
  }

  void unpin_all_nodal_positions()
  {
    unsigned n_node = Fluid_mesh_pt->nnode();
    for (unsigned inod = 0; inod < n_node; inod++)
    {
      Node* nod_pt = Fluid_mesh_pt->node_pt(inod);
      SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);

      solid_node_pt->unpin_position(0);
      solid_node_pt->unpin_position(1);
    }
    std::cout << "Unpin nodal positions " << std::endl;
  }


  void pin_bubble_pressure()
  {
    for (unsigned ibound = First_bubble_boundary_id;
         ibound < First_bubble_boundary_id + 2;
         ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
        nod_pt->pin(0);
      }
    }
    std::cout << "Pin bubble pressure" << std::endl;
  }
  void pin_all_pressure()
  {
    unsigned n_node = Fluid_mesh_pt->nnode();
    for (unsigned inod = 0; inod < n_node; inod++)
    {
      Fluid_mesh_pt->node_pt(inod)->pin(0);
    }
    // Setup equation numbering scheme
    std::cout << "Pin all pressure " << std::endl;
  }
  void unpin_all_pressure()
  {
    unsigned n_node = Fluid_mesh_pt->nnode();
    for (unsigned inod = 0; inod < n_node; inod++)
    {
      Fluid_mesh_pt->node_pt(inod)->unpin(0);
    }
    // Setup equation numbering scheme
    std::cout << "Pin all pressure " << std::endl;
  }

  void pin_curvature_equations()
  {
    //        unsigned nbound=Fluid_mesh_pt->nboundary();
    for (unsigned ibound = First_bubble_boundary_id;
         ibound < First_bubble_boundary_id + 2;
         ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
        nod_pt->pin(1); /// Normal lagrange multiplier
        nod_pt->pin(2); /// Curvature
        nod_pt->pin(3); /// Curvature
        //                nod_pt->pin(4); /// Tangential lagrange multiplier
      }
    }
    std::cout << "Pin curvature equations " << std::endl;
  }

  void pin_tangential_lagrange()
  {
    //        unsigned nbound=Fluid_mesh_pt->nboundary();
    for (unsigned ibound = First_bubble_boundary_id;
         ibound < First_bubble_boundary_id + 2;
         ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
        ////                nod_pt->pin(1); /// Normal lagrange multiplier
        ////                nod_pt->pin(2); /// Curvature
        ////                nod_pt->pin(3); /// Curvature
        nod_pt->pin(4); /// Tangential lagrange multiplier
      }
    }
    std::cout << "Pin tangential lagrange multiplier " << std::endl;
  }

  void unpin_curvature_equations()
  {
    //        unsigned nbound=Fluid_mesh_pt->nboundary();
    for (unsigned ibound = First_bubble_boundary_id;
         ibound < First_bubble_boundary_id + 2;
         ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
        nod_pt->pin(1); /// Normal lagrange multiplier
        nod_pt->unpin(2); /// Curvature
        nod_pt->unpin(3); /// Curvature
        nod_pt->pin(4); /// Tangential lagrange multiplier
      }
    }
    std::cout << "UnPin curvature equations " << std::endl;
  }

  void reset_interface_equations()
  {
    //        unsigned nbound=Fluid_mesh_pt->nboundary();
    for (unsigned ibound = First_bubble_boundary_id;
         ibound < First_bubble_boundary_id + 2;
         ibound++)
    {
      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
        nod_pt->unpin(1); /// Normal lagrange multiplier
        nod_pt->unpin(2); /// Curvature
        nod_pt->unpin(3); /// Curvature
        nod_pt->pin(4); /// Tangential lagrange multiplier
      }
    }
    std::cout << "UnPin curvature equations " << std::endl;
  }

  void redo_equation_numbering()
  {
    cout << "Number of equations: " << this->assign_eqn_numbers() << std::endl;
  }


  //// Operations on parameters

  /// Volume
  void set_V(double new_V)
  {
    Problem_Parameter::Volume = -new_V;
  }
  double get_V()
  {
    return -Problem_Parameter::Volume;
  }
  void increment_V(double dV)
  {
    set_V(get_V() + dV);
  }

  /// Frame speed (if pinned, frame does not move with bubble)
  void set_U(double new_U)
  {
    U_data_pt->set_value(0, new_U);
  }
  double get_U()
  {
    return U_data_pt->value(0);
  }
  void pin_U()
  {
    U_data_pt->pin(0);
  }
  void unpin_U()
  {
    U_data_pt->unpin(0);
  }

  /// Bubble pressure  (volume constraint is not enforced if P is pinned)
  void set_P(double new_P)
  {
    Bubble_pressure_data_pt->set_value(0, new_P);
  }
  double get_P()
  {
    return Bubble_pressure_data_pt->value(0);
  }
  void pin_P()
  {
    Bubble_pressure_data_pt->pin(0);
  }
  void unpin_P()
  {
    Bubble_pressure_data_pt->unpin(0);
  }

  /// Imposed flow rate
  void set_Q(double new_Q)
  {
    Q_inv_data_pt->set_value(0, 1.0 / new_Q);
  }
  double get_Q()
  {
    return 1.0 / Q_inv_data_pt->value(0);
  }
  void increment_Q(double d_Q)
  {
    double new_Q = get_Q() + d_Q;
    set_Q(new_Q);
  }
  double get_Q_inv()
  {
    return Q_inv_data_pt->value(0);
  }

  /// We measure CoM_Y through an integral constraint.
  double get_CoM_Y()
  {
    return CoM_Y_data_pt->value(0);
  }
  void set_CoM_Y(double new_CoM_Y)
  {
    CoM_Y_data_pt->set_value(0, new_CoM_Y);
  }


  /// Channel geometry parameters
  void set_asymmetry(double new_asymmetry)
  {
    Asymmetry_data_pt->set_value(0, new_asymmetry);
  }
  void set_h(double new_h)
  {
    *Problem_Parameter::global_Obstacle_height_pt = new_h;
  }
  void set_w(double new_w)
  {
    *Problem_Parameter::global_Obstacle_width_pt = new_w;
  }
  void set_alpha(double new_alpha)
  {
    *Problem_Parameter::global_alpha_pt = new_alpha;
  }
  double get_asymmetry()
  {
    return Asymmetry_data_pt->value(0);
  }
  double get_h()
  {
    return *Problem_Parameter::global_Obstacle_height_pt;
  }
  double get_w()
  {
    return *Problem_Parameter::global_Obstacle_width_pt;
  }
  double get_maj_rad()
  {
    return Problem_Parameter::Major_Radius;
  }
  double get_min_rad()
  {
    return Problem_Parameter::Minor_Radius;
  }
  double get_xcenter()
  {
    return Problem_Parameter::xcenter;
  }
  double get_ycenter()
  {
    return Problem_Parameter::ycenter;
  }
  unsigned get_circpts()
  {
    return Problem_Parameter::circpts;
  }
  double get_alpha()
  {
    return *Problem_Parameter::global_alpha_pt;
  }
  void increment_h(double dh)
  {
    set_h(get_h() + dh);
  }


  void output_jacobian()
  {
    DoubleVector temp_res;
    CRDoubleMatrix temp_J(this->dof_distribution_pt());
    int n_dof = this->ndof();
    this->get_jacobian(temp_res, temp_J);
    Problem_Parameter::M_file.open("Unstructured_J.dat");
    for (int i = 0; i < n_dof; i++)
    {
      for (int j = 0; j < n_dof; j++)
      {
        Problem_Parameter::M_file << temp_J(i, j) << " ";
      }
      Problem_Parameter::M_file << std::endl;
    }
    Problem_Parameter::M_file.close();
  }

  void set_steady()
  {
    // Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();

    // Loop over them all and make them (temporarily) static
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->make_steady();
    }
  }


private:
  /// \short Create free surface elements
  void create_free_surface_elements();

  /// \short Delete free surface elements
  void delete_free_surface_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = Free_surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Kill surface element
      delete Free_surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Free_surface_mesh_pt->flush_element_and_node_storage();

  } // end of delete_free_surface_elements

  void create_volume_constraint_elements();
  void delete_volume_constraint_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = Volume_constraint_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 1; e < n_element; e++)
    {
      // Kill surface element
      delete Volume_constraint_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Volume_constraint_mesh_pt->flush_element_and_node_storage();
  }

  void create_CoM_X_constraint_elements();
  void delete_CoM_X_constraint_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = CoM_X_constraint_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 1; e < n_element; e++)
    {
      // Kill surface element
      delete CoM_X_constraint_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    CoM_X_constraint_mesh_pt->flush_element_and_node_storage();
  }

  void create_CoM_Y_constraint_elements();
  void delete_CoM_Y_constraint_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = CoM_Y_constraint_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 1; e < n_element; e++)
    {
      // Kill surface element
      delete CoM_Y_constraint_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    CoM_Y_constraint_mesh_pt->flush_element_and_node_storage();
  }


  /// Pointers to mesh of free surface elements
  Mesh* Free_surface_mesh_pt;

  // /// Pointer to mesh containing elements that impose volume constraint
  Mesh* Volume_constraint_mesh_pt;
  Mesh* CoM_X_constraint_mesh_pt;
  Mesh* CoM_Y_constraint_mesh_pt;
  Mesh* Integral_measures_mesh_pt;

  /// Pointer to Fluid_mesh
  RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;

  /// Vector storing pointer to the bubble polygons
  Vector<TriangleMeshPolygon*> Bubble_polygon_pt;

  /// Triangle mesh polygon for outer boundary
  TriangleMeshPolygon* Outer_boundary_polyline_pt;

  Vector<TriangleMeshPolygon*> Obstacle_polygon_pt;

  /// Pointer to a global bubble pressure datum
  Data* Bubble_pressure_data_pt;

  Data* St_data_pt;
  Data* Alpha_data_pt;
  Data* Obstacle_height_data_pt;
  Data* Obstacle_width_data_pt;
  Data* G_data_pt;
  ///---------------------
public:
  Data* U_data_pt;
  Data* Q_inv_data_pt;
  ///---------------------
  Data* Asymmetry_data_pt;
  Data* CoM_Y_data_pt;

  Data* Integral_measures_data_pt;

  Mesh* Inflow_mesh_pt;
  Mesh* Outflow_mesh_pt;
  //// /// Pointer to element that imposes volume constraint for bubble
  VolumeConstraintElement* Vol_constraint_el_pt;
  VolumeConstraintElement* CoM_X_constraint_el_pt;
  //    VolumeConstraintElement* CoM_Y_constraint_el_pt;
  SelfReferentialVolumeConstraintElement* CoM_Y_constraint_el_pt;
  SelfReferentialVolumeConstraintElement* Integral_measures_el_pt;

  /// Enumeration of mesh boundaries
  enum
  {
    Inflow_boundary_id = 0,
    Upper_wall_boundary_id = 1,
    Outflow_boundary_id = 2,
    Bottom_wall_boundary_id = 3,
    First_bubble_boundary_id = 4,
    Second_bubble_boundary_id = 5,
  };


  ///-------------------------------------------------------------------------------

  void set_steady_constraints()
  {
    std::cout << "SET STEADY CONSTRAINTS !!!!!!!!!!" << std::endl;
    unpin_all_pressure();
    unpin_all_nodal_positions();
    // unpin_lagrange_multipliers();
    unpin_curvature_equations();
    // pin_corner_curvature_equations();
    // pin_tangential_lagrange_multipliers();
    // G_at_upper_outflow_data_pt->pin(0);
    // G_at_upper_outflow_data_pt->set_value(0,0);
    //        remove_corner_curvature();
    // Problem_Parameter::steady_version_pt = true;

    //        set_tip_node_pt();

    /// Now, set wall conditions.
    complete_problem_setup();

    redo_equation_numbering();
  }


  void set_steady_constraints_free_Q()
  {
    std::cout << "========================" << std::endl;
    std::cout << "Set steady constraints free Q" << std::endl;
    std::cout << "G is pinned. Q varies to meet position constraint"
              << std::endl;
    std::cout << "========================" << std::endl;

    set_steady_constraints();
    pin_G();
    unpin_Q();
    pin_U();

    redo_equation_numbering();
  }


  void pin_G()
  {
    G_data_pt->pin(0);
  }
  void pin_Q()
  {
    Q_inv_data_pt->pin(0);
  }

  void unpin_G()
  {
    G_data_pt->unpin(0);
  }
  void unpin_Q()
  {
    Q_inv_data_pt->unpin(0);
  }

  bool get_jacobian_sign_change()
  {
    return First_jacobian_sign_change;
  }
  ///-------------------------------------------------------------------------------
}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
BubbleInChannelProblem<ELEMENT>::BubbleInChannelProblem()
{
  ///=============================
  Bifurcation_detection = true;
  Bisect_to_find_bifurcation = false;
  Desired_newton_iterations_ds = 100;
  ///=============================

  // Allocate the timestepper -- this constructs the Problem's
  // time object with a sufficient amount of storage to store the
  // previous timsteps.
  bool adaptive_timestepping = true;
  add_time_stepper_pt(new BDF<2>(adaptive_timestepping));
  double initial_height = 0.00;
  double initial_width = 0.25;
  double initial_alpha = 10;
  double initial_St = 1;
  double initial_G = 1;
  double initial_Q_inv = 20;
  double initial_Frame_speed = 0;
  double initial_asymmetry = 0.00;
  double initial_CoM_Y = 0.0; /// initially 0
  double initial_bubble_pressure = 0;

  /// Or might vary Q.
  Q_inv_data_pt = new Data(1);
  this->add_global_data(Q_inv_data_pt);
  Q_inv_data_pt->set_value(0, initial_Q_inv);
  Problem_Parameter::global_Q_inv_pt = Q_inv_data_pt->value_pt(0);
  Q_inv_data_pt->pin(0);

  /// Obstacle height
  Obstacle_height_data_pt = new Data(1);
  this->add_global_data(Obstacle_height_data_pt);
  Obstacle_height_data_pt->set_value(0, initial_height);
  Problem_Parameter::global_Obstacle_height_pt =
    Obstacle_height_data_pt->value_pt(0);
  Obstacle_height_data_pt->pin(0);

  /// Obstacle width
  Obstacle_width_data_pt = new Data(1);
  this->add_global_data(Obstacle_width_data_pt);
  Obstacle_width_data_pt->set_value(0, initial_width);
  Problem_Parameter::global_Obstacle_width_pt =
    Obstacle_width_data_pt->value_pt(0);
  Obstacle_width_data_pt->pin(0);

  /// Alpha (should be greater than 1)
  Alpha_data_pt = new Data(1);
  this->add_global_data(Alpha_data_pt);
  Alpha_data_pt->set_value(0, initial_alpha);
  Problem_Parameter::global_alpha_pt = Alpha_data_pt->value_pt(0);
  Alpha_data_pt->pin(0);

  /// Strouhal number
  St_data_pt = new Data(1);
  this->add_global_data(St_data_pt);
  St_data_pt->set_value(0, initial_St);
  Problem_Parameter::global_St_pt = St_data_pt->value_pt(0);
  St_data_pt->pin(0);

  CoM_Y_data_pt = new Data(1);
  this->add_global_data(CoM_Y_data_pt);
  CoM_Y_data_pt->set_value(0, initial_CoM_Y);
  Problem_Parameter::global_CoM_Y_pt = CoM_Y_data_pt->value_pt(0);
  CoM_Y_data_pt->unpin(0);

  /// Position can be constrained by varying G: a pressure gradient
  G_data_pt = new Data(1);
  this->add_global_data(G_data_pt);
  G_data_pt->set_value(0, initial_G);
  Problem_Parameter::global_G_pt = G_data_pt->value_pt(0);
  G_data_pt->pin(0);

  /// Asymmetry parameter
  Asymmetry_data_pt = new Data(1);
  this->add_global_data(Asymmetry_data_pt);
  Asymmetry_data_pt->set_value(0, initial_asymmetry);
  Problem_Parameter::global_Asymmetry_pt = Asymmetry_data_pt->value_pt(0);
  Asymmetry_data_pt->pin(0);


  /// Bubble pressure
  Bubble_pressure_data_pt = new Data(1);
  this->add_global_data(Bubble_pressure_data_pt);
  Bubble_pressure_data_pt->set_value(0, initial_bubble_pressure);
  Problem_Parameter::global_Bubble_pressure_pt =
    Bubble_pressure_data_pt->value_pt(0);
  Bubble_pressure_data_pt->unpin(0);


  U_data_pt = new Data(1);
  this->add_global_data(U_data_pt);
  U_data_pt->set_value(0, initial_Frame_speed);
  Problem_Parameter::global_Frame_speed_pt = U_data_pt->value_pt(0);
  U_data_pt->unpin(0);

  Integral_measures_data_pt = new Data(Problem_Parameter::n_integral_measures);
  this->add_global_data(Integral_measures_data_pt);
  for (unsigned i = 0; i < Problem_Parameter::n_integral_measures; i++)
  {
    Integral_measures_data_pt->set_value(i, 0.0);
    Integral_measures_data_pt->unpin(i);
  }


  unsigned index_of_traded_pressure = 0;
  Vol_constraint_el_pt = new VolumeConstraintElement(&Problem_Parameter::Volume,
                                                     Bubble_pressure_data_pt,
                                                     index_of_traded_pressure);

  CoM_X_constraint_el_pt = new VolumeConstraintElement(
    &Problem_Parameter::Centre_of_mass, U_data_pt, index_of_traded_pressure);

  CoM_Y_constraint_el_pt = new SelfReferentialVolumeConstraintElement(
    Problem_Parameter::global_CoM_Y_pt,
    CoM_Y_data_pt,
    index_of_traded_pressure);

  Integral_measures_el_pt = new SelfReferentialVolumeConstraintElement(
    Problem_Parameter::global_CoM_Y_pt,
    Integral_measures_data_pt,
    index_of_traded_pressure);


  // Build the boundary segments for outer boundary, consisting of
  //--------------------------------------------------------------
  // four separate polylines
  //------------------------
  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

  // Each polyline only has two vertices -- provide storage for their
  // coordinates
  Vector<Vector<double>> vertex_coord(2);
  for (unsigned i = 0; i < 2; i++)
  {
    vertex_coord[i].resize(2);
  }

  // First polyline: Inflow
  vertex_coord[0][0] = -Problem_Parameter::Length;
  vertex_coord[0][1] = -1.0;
  vertex_coord[1][0] = -Problem_Parameter::Length;
  vertex_coord[1][1] = 1.0;

  // Build the 1st boundary polyline
  boundary_polyline_pt[0] =
    new TriangleMeshPolyLine(vertex_coord, Inflow_boundary_id);

  // Second boundary polyline: Upper wall
  vertex_coord[0][0] = -Problem_Parameter::Length;
  vertex_coord[0][1] = 1.0;
  vertex_coord[1][0] = Problem_Parameter::Length;
  vertex_coord[1][1] = 1.0;

  // Build the 2nd boundary polyline
  boundary_polyline_pt[1] =
    new TriangleMeshPolyLine(vertex_coord, Upper_wall_boundary_id);

  // Third boundary polyline: Outflow
  vertex_coord[0][0] = Problem_Parameter::Length;
  vertex_coord[0][1] = 1.0;
  vertex_coord[1][0] = Problem_Parameter::Length;
  vertex_coord[1][1] = -1.0;

  // Build the 3rd boundary polyline
  boundary_polyline_pt[2] =
    new TriangleMeshPolyLine(vertex_coord, Outflow_boundary_id);

  // Fourth boundary polyline: Bottom wall
  vertex_coord[0][0] = Problem_Parameter::Length;
  vertex_coord[0][1] = -1.0;
  vertex_coord[1][0] = -Problem_Parameter::Length;
  vertex_coord[1][1] = -1.0;

  // Build the 4th boundary polyline
  boundary_polyline_pt[3] =
    new TriangleMeshPolyLine(vertex_coord, Bottom_wall_boundary_id);


  // Create the triangle mesh polygon for outer boundary
  Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);


  // Now define initial shape of bubble(s) with polygon
  //---------------------------------------------------

  // We have one bubble
  Bubble_polygon_pt.resize(1);

  for (unsigned i_circle = 0; i_circle < 1; i_circle++)

  {
    // Place it smack in the middle of the channel
    double x_center = Problem_Parameter::xcenter;
    double y_center = Problem_Parameter::ycenter;
    double major_radius = Problem_Parameter::Major_Radius;
    double minor_radius = Problem_Parameter::Minor_Radius;
    Ellipse* bubble_pt = new Ellipse(major_radius, minor_radius);

    // Intrinsic coordinate along GeomObject defining the bubble
    Vector<double> zeta(1);

    // Position vector to GeomObject defining the bubble
    Vector<double> coord(2);

    // Number of points defining bubble
    unsigned npoints = Problem_Parameter::circpts; // 16;
    double unit_zeta = MathematicalConstants::Pi / double(npoints - 1);

    // This bubble is bounded by two distinct boundaries, each
    // represented by its own polyline
    Vector<TriangleMeshCurveSection*> bubble_polyline_pt(2);

    // Vertex coordinates
    Vector<Vector<double>> bubble_vertex(npoints);

    // Create points on boundary
    for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
    {
      // Resize the vector
      bubble_vertex[ipoint].resize(2);

      // Get the coordinates
      zeta[0] = unit_zeta * double(ipoint);
      bubble_pt->position(zeta, coord);

      // Shift
      bubble_vertex[ipoint][0] = coord[0] + x_center;
      bubble_vertex[ipoint][1] = coord[1] + y_center;
    }

    // Build the 1st bubble polyline
    bubble_polyline_pt[0] = new TriangleMeshPolyLine(
      bubble_vertex, First_bubble_boundary_id + 2 * i_circle);

    // Second boundary of bubble
    for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
    {
      // Resize the vector
      bubble_vertex[ipoint].resize(2);

      // Get the coordinates
      zeta[0] = (unit_zeta * double(ipoint)) + MathematicalConstants::Pi;
      bubble_pt->position(zeta, coord);

      // Shift
      bubble_vertex[ipoint][0] = coord[0] + x_center;
      bubble_vertex[ipoint][1] = coord[1] + y_center;
    }

    // Build the 2nd bubble polyline
    bubble_polyline_pt[1] = new TriangleMeshPolyLine(
      bubble_vertex, Second_bubble_boundary_id + 2 * i_circle);


    // Define coordinates of a point inside the bubble
    Vector<double> bubble_center(2);
    bubble_center[0] = x_center;
    bubble_center[1] = y_center;

    bubble_polyline_pt[0]->set_maximum_length(0.02);
    bubble_polyline_pt[1]->set_maximum_length(0.02);

    // Create closed polygon from two polylines
    Bubble_polygon_pt[i_circle] =
      new TriangleMeshPolygon(bubble_polyline_pt, bubble_center);


    Bubble_polygon_pt[i_circle]->set_polyline_refinement_tolerance(0.05);
    Bubble_polygon_pt[i_circle]->set_polyline_unrefinement_tolerance(0.01);
    std::cout << "Bubble " << i_circle << std::endl;
  }


  // Now build the mesh, based on the boundaries specified by
  //---------------------------------------------------------
  // polygons just created
  //----------------------

  // Convert to "closed curve" objects
  TriangleMeshClosedCurve* outer_closed_curve_pt = Outer_boundary_polyline_pt;
  unsigned nb = Bubble_polygon_pt.size();
  Vector<TriangleMeshClosedCurve*> bubble_closed_curve_pt(nb);
  for (unsigned i = 0; i < nb; i++)
  {
    bubble_closed_curve_pt[i] = Bubble_polygon_pt[i];
  }

  // Target area for initial mesh
  double uniform_element_area = 0.1;

  // Use the TriangleMeshParameters object for gathering all
  // the necessary arguments for the TriangleMesh object
  TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

  // Define the holes on the boundary
  triangle_mesh_parameters.internal_closed_curve_pt() = bubble_closed_curve_pt;

  // Define the maximum element areas
  triangle_mesh_parameters.element_area() = uniform_element_area;

  // Create the mesh
  Fluid_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
    triangle_mesh_parameters, this->time_stepper_pt());

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Fluid_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  Fluid_mesh_pt->max_permitted_error() = 2e-5; /// 5e-2;
  Fluid_mesh_pt->min_permitted_error() = 5e-6; /// 1e-4;
  Fluid_mesh_pt->max_element_size() = 0.05; /// 0.2;
  Fluid_mesh_pt->min_element_size() = 1e-6; /// 4e-5;

  // Set boundary condition and complete the build of all elements
  complete_problem_setup();

  // Construct the mesh of free surface elements
  Free_surface_mesh_pt = new Mesh;
  create_free_surface_elements();

  Volume_constraint_mesh_pt = new Mesh;
  create_volume_constraint_elements();


  CoM_X_constraint_mesh_pt = new Mesh;
  create_CoM_X_constraint_elements();


  CoM_Y_constraint_mesh_pt = new Mesh;
  create_CoM_Y_constraint_elements();

  Inflow_mesh_pt = new Mesh;
  create_inflow_elements();

  Integral_measures_mesh_pt = new Mesh;
  Integral_measures_mesh_pt->add_element_pt(Integral_measures_el_pt);

  //    Outflow_mesh_pt = new Mesh;
  //    create_outflow_elements();

  // Add Fluid_mesh_pt sub meshes
  this->add_sub_mesh(Fluid_mesh_pt);

  // Add Free_surface sub meshes
  this->add_sub_mesh(this->Free_surface_mesh_pt);
  this->add_sub_mesh(this->Volume_constraint_mesh_pt);
  this->add_sub_mesh(this->CoM_X_constraint_mesh_pt);
  this->add_sub_mesh(this->CoM_Y_constraint_mesh_pt);
  this->add_sub_mesh(this->Integral_measures_mesh_pt);

  this->add_sub_mesh(this->Inflow_mesh_pt);
  //    this->add_sub_mesh(this->Outflow_mesh_pt);

  // Build global mesh
  this->build_global_mesh();

  pin_tangential_lagrange();

  Use_finite_differences_for_continuation_derivatives = true;

  Max_residuals = 300000000000;
  Max_newton_iterations = 300;
  // Use_continuation_timestepper=true;
  Desired_proportion_of_arc_length = 0.00001;
  Desired_newton_iterations_ds = 300;

  // Setup equation numbering scheme
  cout << "Number of equations: " << this->assign_eqn_numbers() << std::endl;
  //  linear_solver_pt()=new FD_LU;
} // end_of_constructor


//============start_of_create_free_surface_elements======================
/// Create elements that impose the kinematic and dynamic bcs
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_free_surface_elements()
{
  // Volume constraint element stores the Data item that stores
  // the bubble pressure that is adjusted/tradto allow for
  // volume conservation. Which value is the pressure stored in?
  //// unsigned p_traded_index=Vol_constraint_el_pt->index_of_traded_pressure();

  // Loop over the free surface boundaries
  //    unsigned nb=Fluid_mesh_pt->nboundary();
  for (unsigned b = First_bubble_boundary_id; b < Second_bubble_boundary_id + 1;
       b++)
  {
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

      HeleShawInterfaceElement<ELEMENT>* el_pt =
        new HeleShawInterfaceElement<ELEMENT>(bulk_elem_pt, face_index);


      // Add it to the mesh
      Free_surface_mesh_pt->add_element_pt(el_pt);

      // Add the appropriate boundary number
      el_pt->set_boundary_number_in_bulk_mesh(b);

      //                /// To calculate constants for surface tension
      el_pt->q_inv_pt() = Q_inv_data_pt->value_pt(0);

      // Set Strouhal number
      el_pt->st_pt() = St_data_pt->value_pt(0);

      // Set alpha
      el_pt->alpha_pt() = Alpha_data_pt->value_pt(0);

      el_pt->upper_wall_fct_pt() = Problem_Parameter::upper_wall_fct;

      //            el_pt->upper_wall_flux_fct_pt() =
      //            Problem_Parameter::upper_wall_flux_fct;

      el_pt->wall_speed_fct_pt() = Problem_Parameter::frame_speed_fct;

      el_pt->bubble_pressure_fct_pt() =
        Problem_Parameter::bubble_pressure_function;
      // jack
      el_pt->add_external_data(Q_inv_data_pt, true);
      el_pt->add_external_data(U_data_pt, true);
      el_pt->add_external_data(Bubble_pressure_data_pt, true);
      el_pt->add_external_data(Alpha_data_pt, true);

      unsigned measure_data_number =
        el_pt->add_external_data(Integral_measures_data_pt, true);
      el_pt->integral_measure_data_number = measure_data_number;

      //            el_pt->use_minisolve_state_pt =
      //            &Problem_Parameter::Minisolve_state;
    }
  }
}
// end of create_free_surface_elements

//============start_of_create_volume_constraint_elements=================
/// Create elements that impose volume constraint on the bubble
//=======================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_volume_constraint_elements()
{
  // Add volume constraint element to the mesh
  Volume_constraint_mesh_pt->add_element_pt(Vol_constraint_el_pt);

  // Loop over the free surface boundaries
  // unsigned nb=Fluid_mesh_pt->nboundary();
  for (unsigned b = First_bubble_boundary_id; b < Second_bubble_boundary_id + 1;
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
        new HeleShawVolumeConstraintElement<ELEMENT>(bulk_elem_pt, face_index);

      // Set the "master" volume constraint element
      el_pt->set_volume_constraint_element(Vol_constraint_el_pt);

      // Add it to the mesh
      Volume_constraint_mesh_pt->add_element_pt(el_pt);
    }
  }
}
// end of create_volume_constraint_elements


//============start_of_create_volume_constraint_elements=================
/// Create elements that impose volume constraint on the bubble
//=======================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_CoM_X_constraint_elements()
{
  // Add volume constraint element to the mesh
  CoM_X_constraint_mesh_pt->add_element_pt(CoM_X_constraint_el_pt);


  for (unsigned b = First_bubble_boundary_id; b < Second_bubble_boundary_id + 1;
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
      HeleShawCoMXConstraintElement<ELEMENT>* el_pt =
        new HeleShawCoMXConstraintElement<ELEMENT>(bulk_elem_pt, face_index);

      // Set the "master" volume constraint element
      el_pt->set_volume_constraint_element(CoM_X_constraint_el_pt);

      // Add it to the mesh
      CoM_X_constraint_mesh_pt->add_element_pt(el_pt);
    }
  }
}
// end of create_volume_constraint_elements

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_CoM_Y_constraint_elements()
{
  // Add volume constraint element to the mesh
  CoM_Y_constraint_mesh_pt->add_element_pt(CoM_Y_constraint_el_pt);


  for (unsigned b = First_bubble_boundary_id; b < Second_bubble_boundary_id + 1;
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
      HeleShawCoMYConstraintElement<ELEMENT>* el_pt =
        new HeleShawCoMYConstraintElement<ELEMENT>(bulk_elem_pt, face_index);

      // Set the "master" volume constraint element
      el_pt->set_volume_constraint_element(CoM_Y_constraint_el_pt);

      // Add it to the mesh
      CoM_Y_constraint_mesh_pt->add_element_pt(el_pt);
    }
  }
}

//==start_of_complete_problem_setup=======================================
/// Set boundary conditions and complete the build of all elements
//========================================================================

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::complete_problem_setup()
{
  // Map to record if a given boundary is on a bubble or not
  map<unsigned, bool> is_on_bubble_bound;

  // Loop over the bubbles
  // unsigned nbubble=Bubble_polygon_pt.size();
  unsigned nbubble = 1;
  for (unsigned ibubble = 0; ibubble < nbubble; ibubble++)
  {
    // Get the vector all boundary IDs associated with the polylines that
    // make up the closed polygon
    Vector<unsigned> bubble_bound_id =
      this->Bubble_polygon_pt[ibubble]->polygon_boundary_id();

    // Get the number of boundary
    unsigned nbound = bubble_bound_id.size();

    // Fill in the map
    for (unsigned ibound = 0; ibound < nbound; ibound++)
    {
      // This boundary...
      unsigned bound_id = bubble_bound_id[ibound];

      // ...is on the bubble
      is_on_bubble_bound[bound_id] = true;
      std::cout << bound_id << " on bubble" << std::endl;
    }
  } // points on bubble boundary located

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
      Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);

      // Pin pseudo-solid positions apart from bubble boundary which
      // we allow to move
      SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
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

      if (ibound == Outflow_boundary_id)
      {
        nod_pt->pin(0);
        nod_pt->set_value(0, 0.0);
        //                std::cout << "Pin pressure to zero on outflow boundary
        //                " << std::endl;
      }
    }
  } // end loop over boundaries

  // Complete the build of all elements so they are fully functional
  // Remember that adaptation for triangle meshes involves a complete
  // regneration of the mesh (rather than splitting as in tree-based
  // meshes where such parameters can be passed down from the father
  // element!)
  unsigned n_element = Fluid_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

    // Set the constitutive law for pseudo-elastic mesh deformation
    el_pt->constitutive_law_pt() = Problem_Parameter::Constitutive_law_pt;
    el_pt->upper_wall_fct_pt() = Problem_Parameter::upper_wall_fct;
  }

  //    cout <<"Number of equations: " << this->assign_eqn_numbers() <<
  //    std::endl;
} // end of complete_problem_setup

//========================================================================
/// Compute error estimates and assign to elements for plotting
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::compute_error_estimate(double& max_err,
                                                             double& min_err)
{
  std::cout << "Compute error estimate" << std::endl;
  // Get error estimator
  ErrorEstimator* err_est_pt = Fluid_mesh_pt->spatial_error_estimator_pt();

  // Get/output error estimates
  unsigned nel = Fluid_mesh_pt->nelement();
  Vector<double> elemental_error(nel);

  // We need a dynamic cast, get_element_errors from the Fluid_mesh_pt
  // Dynamic cast is used because get_element_errors require a Mesh* ans
  // not a SolidMesh*
  Mesh* fluid_mesh_pt = dynamic_cast<Mesh*>(Fluid_mesh_pt);
  err_est_pt->get_element_errors(fluid_mesh_pt, elemental_error);

  // Set errors for post-processing and find extrema
  max_err = 0.0;
  min_err = 1e6;
  for (unsigned e = 0; e < nel; e++)
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
void BubbleInChannelProblem<ELEMENT>::jack_solve(double Q,
                                                 bool jack_arc,
                                                 unsigned n_iter,
                                                 double sign)
{
  if (jack_arc == true)
  {
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "PSEUDO-ARCLENGTH CONTINUATION STARTING=====" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
  }
  else
  {
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "NORMAL-ARCLENGTH CONTINUATION STARTING=====" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
  }

  reset_arc_length_parameters();
  Desired_newton_iterations_ds = 3;
  Desired_proportion_of_arc_length = 0.1;
  double new_Q = Q;
  double new_Q_inv = 1 / Q;
  double ds;
  if (jack_arc == true)
  {
    ds = -1; //-10
  }
  else
  {
    ds = -1; //-10
  }
  ds = -1;
  for (unsigned m = 0; m < n_iter; m++)
  {
    Q_inv_data_pt->pin(0);

    //      assign_eqn_numbers();
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "CONTINUATION STEP " << m << std::endl;
    std::cout << "U_b is " << get_U() << std::endl;
    std::cout << "Ca is " << get_Q() * get_U() << std::endl;
    std::cout << "ds is " << ds << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    if (m % 5 == 4)
    {
      if (jack_arc == true)
      {
        ds = arc_length_step_solve(Q_inv_data_pt, 0, ds);
        reset_arc_length_parameters();
        steady_newton_solve(0);
        steady_newton_solve(1);
        Q_inv_data_pt->pin(0);
      }
      else
      {
        steady_newton_solve(0);
        steady_newton_solve(1);
      }
    }
    else
    {
      if (jack_arc == true)
      {
        ds = arc_length_step_solve(Q_inv_data_pt, 0, ds);
      }
      else
      {
        steady_newton_solve(0);
        new_Q_inv += ds;
        new_Q = 1 / new_Q_inv;
        set_Q(new_Q);
      }
    }
    reset_lagrangian_coordinates();
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "COM_Y is " << get_CoM_Y() << std::endl;
    std::cout << "U_b is " << get_U() << std::endl;
    std::cout << "Ca is " << get_U() * get_Q() << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    doc_solution();
  }
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::dump_it(ofstream& dump_file)
{
  dump(dump_file);
}
