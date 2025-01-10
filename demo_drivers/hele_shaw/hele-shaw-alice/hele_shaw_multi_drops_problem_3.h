
//==start_of_problem_class============================================
/// Problem class to simulate inviscid bubble propagating along 2D channel
//====================================================================
template<class ELEMENT>
class BubbleInChannelProblem : public Problem
{
public:
  /// Constructor
  BubbleInChannelProblem();

  /// Destructor
  ~BubbleInChannelProblem()
  {
    delete this->time_stepper_pt(0);

    unsigned n = Outer_boundary_polyline_pt->npolyline();
    for (unsigned j = 0; j < n; j++)
    {
      delete Outer_boundary_polyline_pt->polyline_pt(j);
    }
    delete Outer_boundary_polyline_pt;

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


    delete_free_surface_elements();
    delete Free_surface_mesh_pt;

    delete_inflow_elements();
    delete Inflow_mesh_pt;

    delete_CoM_X_constraint_elements();
    delete CoM_X_constraint_mesh_pt;

    delete_CoM_Y_constraint_elements();
    delete CoM_Y_constraint_mesh_pt;

    delete_volume_constraint_elements();
    delete Volume_constraint_mesh_pt;

    delete Fluid_mesh_pt->spatial_error_estimator_pt();

    delete Fluid_mesh_pt;
  }


  void actions_before_adapt()
  {
    // if (Constraints_are_steady)
    //{std::cout << "In actions before adapt, with steady constraints"<<
    //std::endl;} else {std::cout << "In actions before adapt, with unsteady
    //constraints"<< std::endl;}

    delete_free_surface_elements();
    delete_inflow_elements();
    delete_CoM_X_constraint_elements();
    delete_CoM_Y_constraint_elements();
    //		delete Free_surface_mesh_pt;
    //		delete Inflow_mesh_pt;
    //		delete CoM_X_constraint_mesh_pt;

    if (Constraints_are_steady)
    {
      delete_volume_constraint_elements();
      //			delete Volume_constraint_mesh_pt;
    }

    //		this->flush_sub_meshes();
    this->rebuild_global_mesh();
    //		this->add_sub_mesh(Fluid_mesh_pt);
    //		this->rebuild_global_mesh();
    // redo_equation_numbering();
  }


  void actions_after_adapt()
  {
    if (Constraints_are_steady)
    {
      std::cout << "In actions after adapt, with steady constraints"
                << std::endl;
    }
    else
    {
      std::cout << "In actions after adapt, with unsteady constraints"
                << std::endl;
    }

    complete_problem_setup();

    //		Free_surface_mesh_pt=new Mesh;
    create_free_surface_elements();
    //		Inflow_mesh_pt = new Mesh;
    create_inflow_elements();
    if (Constraints_are_steady)
    {
      //		Volume_constraint_mesh_pt = new Mesh;
      create_volume_constraint_elements();
    }
    //		CoM_X_constraint_mesh_pt = new Mesh;
    create_CoM_X_constraint_elements();
    create_CoM_Y_constraint_elements();

    //		this->add_sub_mesh(Fluid_mesh_pt);
    //		this->add_sub_mesh(this->Free_surface_mesh_pt);
    //		this->add_sub_mesh(this->Inflow_mesh_pt);
    //		this->add_sub_mesh(this->Volume_constraint_mesh_pt);
    //		this->add_sub_mesh(this->CoM_X_constraint_mesh_pt);

    this->rebuild_global_mesh();


    redo_equation_numbering();
  }


  double global_temporal_error_norm()
  {
    //	std::cout << "Calling global error norm" << std::endl;
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
  void actions_after_newton_solve()
  {
    // check_topology();
  }

  void actions_after_newton_step()
  {
    if (Problem_Parameter::doc_every_step == true)
    {
      doc_solution();
    }
  }

  /// Update the problem specs before solve
  void actions_before_newton_solve() {}

  // void actions_before_newton_convergence_check(){doc_solution();}

  void complete_problem_setup();

  void doc_solution(const std::string& comment = "");

  void compute_error_estimate(double& max_err, double& min_err);

  void reset_lagrangian_coordinates()
  {
    Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
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
    //	cout << n_inflow_element << " Poisson Flux Elements created on main
    //inflow boundary" << endl; 	std::cout << "G is "<< G_data_pt->value(0)
    //<<std::endl;
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

  void create_CoM_X_constraint_elements();
  void delete_CoM_X_constraint_elements()
  {
    unsigned n_element = CoM_X_constraint_mesh_pt->nelement();
    for (unsigned e = 1; e < n_element; e++)
    {
      delete CoM_X_constraint_mesh_pt->element_pt(e);
    }
    CoM_X_constraint_mesh_pt->flush_element_and_node_storage();
  }


  void redo_equation_numbering()
  {
    this->assign_eqn_numbers();
    //	cout <<"Number of equations: " << this->assign_eqn_numbers() <<
    //std::endl;
  }


  //// Operations on parameters

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

  void set_G(double new_G)
  {
    G_data_pt->set_value(0, new_G);
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
  void set_drop_viscosity(double new_drop_viscosity)
  {
    *Problem_Parameter::global_drop_viscosity_pt = new_drop_viscosity;
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
  double get_alpha()
  {
    return *Problem_Parameter::global_alpha_pt;
  }
  double get_drop_viscosity()
  {
    return *Problem_Parameter::global_drop_viscosity_pt;
  }

  void increment_h(double dh)
  {
    set_h(get_h() + dh);
  }

  double get_time()
  {
    return this->time_pt()->time();
  }

  void set_steady()
  {
    // std::cout << "Running set steady!!! " <<std::endl;
    // Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();

    // Loop over them all and make them (temporarily) static
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->make_steady();
    }
  }

  void write_output_file();

  void solve_for_eigenproblem();


  void assess_bubble_volumes(bool reset_bubble_volume);


  void snap_bubble_to_ellipse(
    unsigned i_bubble, double x_center, double y_center, double a, double b)
  {
    /// Travel over all nodes on the boundary of bubble i_bubble, and move
    /// radially outwards onto an ellipse.

    /// Each bubble is constructed from two boundaries.
    for (unsigned i_b = 0; i_b < 2; i_b++)
    {
      unsigned ibound = 4 + 2 * (i_bubble) + i_b;

      unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);

        SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);

        double current_x = solid_node_pt->position(0);
        double current_y = solid_node_pt->position(1);

        double current_distance =
          std::sqrt((current_x - x_center) * (current_x - x_center) / (a * a) +
                    (current_y - y_center) * (current_y - y_center) / (b * b));

        double new_x = x_center + (current_x - x_center) / current_distance;
        double new_y = y_center + (current_y - y_center) / current_distance;

        solid_node_pt->x(0) = new_x;
        solid_node_pt->x(1) = new_y;
        //	std::cout << solid_node_pt->value(0)<< " " <<
        //solid_node_pt->position(0) << std::endl;
      } /// End loop over nodes


    } /// End loop over the two boundaries on this bubble
  }

  void output_residuals()
  {
    DoubleVector temp_res;
    int n_dof = this->ndof();
    CRDoubleMatrix temp_J(this->dof_distribution_pt());
    this->get_jacobian(temp_res, temp_J);
    Problem_Parameter::M_file.open("Residuals.dat");
    for (int i = 0; i < n_dof; i++)
    {
      Problem_Parameter::M_file << temp_res[i] << std::endl;
    }
    Problem_Parameter::M_file.close();
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


  void output_eigenmatrices()
  {
    std::cout << "OUTPUT EIGENMATRICES" << std::endl;
    CRDoubleMatrix temp_M(this->dof_distribution_pt());
    CRDoubleMatrix temp_J(this->dof_distribution_pt());

    Problem::get_eigenproblem_matrices(temp_M, temp_J, 0.0);

    std::cout << "Now write to file" << std::endl;

    int n_dof = this->ndof();
    /// Write the original A and M to file.
    Problem_Parameter::M_file.open("Unstructured_M.dat");
    for (int i = 0; i < n_dof; i++)
    {
      for (int j = 0; j < n_dof; j++)
      {
        Problem_Parameter::M_file << temp_M(i, j) << " ";
      }
      Problem_Parameter::M_file << std::endl;
    }
    Problem_Parameter::M_file.close();


    /// Write the original A and M to file.
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

  void set_steady_constraints()
  {
    if (Constraints_are_steady == false)
    {
      // std::cout << "Currently created as unsteady, so create elements" <<
      // std::endl;
      create_volume_constraint_elements();
    }
    Constraints_are_steady = true;
    // std::cout << "Rebuild global mesh" << std::endl;
    this->rebuild_global_mesh();
    redo_equation_numbering();
    // cout <<"Number of equations: " << this->assign_eqn_numbers() <<
    // std::endl;
  }

  void set_unsteady_constraints()
  {
    if (Constraints_are_steady == true)
    {
      std::cout << "Delete the volume constraint elements... " << std::endl;
      delete_volume_constraint_elements();
    }
    Constraints_are_steady = false;
    this->rebuild_global_mesh();
    redo_equation_numbering();
    std::cout << "Have now set unsteady constraints" << std::endl;
    // cout <<"Number of equations: " << this->assign_eqn_numbers() <<
    // std::endl;
  }

  bool Constraints_are_steady;

  void output_ordered_polygon(ofstream& Boundary_output_file)
  {
    Problem_Parameter::Ordered_bubbles.resize(Problem_Parameter::N_Bubble);
    unsigned N_Bubble_initial = Problem_Parameter::N_Bubble;

    /// Get an ordered representation of each bubble
    for (unsigned bubble = 0; bubble < N_Bubble_initial; bubble++)
    {
      unsigned num_nodes = Fluid_mesh_pt->nboundary_node(4 + 2 * (bubble));
      num_nodes += Fluid_mesh_pt->nboundary_node(4 + 2 * (bubble) + 1);

      Problem_Parameter::Ordered_bubbles[bubble].resize(num_nodes);

      bool travel_forward = true;
      unsigned bubble_node = 0;

      for (unsigned i_b = 0; i_b < 2; i_b++)
      {
        if (i_b == 0)
        {
          travel_forward = true;
        }
        else
        {
          travel_forward = false;
        }
        unsigned b = 4 + 2 * (bubble) + i_b;

        std::set<Node*> boundary_nodes_pt;
        const unsigned n_boundary_ele = Fluid_mesh_pt->nboundary_element(b);
        for (unsigned e = 0; e < n_boundary_ele; e++)
        {
          // Get the boundary bulk element
          FiniteElement* bulk_ele_pt = Fluid_mesh_pt->boundary_element_pt(b, e);
          // Get the face index
          int face_index = Fluid_mesh_pt->face_index_at_boundary(b, e);
          // Create the face element
          FiniteElement* face_ele_pt =
            new DummyFaceElement<ELEMENT>(bulk_ele_pt, face_index);

          // Get the number of nodes on the face element
          const unsigned n_nodes = face_ele_pt->nnode();
          for (unsigned i = 0; i < n_nodes; i++)
          {
            // Get the nodes in the face elements
            Node* tmp_node_pt = face_ele_pt->node_pt(i);
            // Add the nodes to the set of boundary nodes
            boundary_nodes_pt.insert(tmp_node_pt);
          } // for (i < n_nodes)

          // Free the memory allocated for the face element
          delete face_ele_pt;
          face_ele_pt = 0;
        } // for (e < n_boundary_ele)

        // Set to store the boundary nodes in order
        std::set<Vector<double>> set_node_coord;
        // Loop over the nodes on the boundary and store them in the set
        for (std::set<Node*>::iterator it = boundary_nodes_pt.begin();
             it != boundary_nodes_pt.end();
             it++)
        {
          Node* inode_pt = (*it);

          // Get the node coordinates
          const unsigned n_dim = inode_pt->ndim();
          Vector<double> node_coord(n_dim + 1);

          // Get the boundary coordinate
          Vector<double> zeta(1);
          inode_pt->get_coordinates_on_boundary(b, zeta);
          node_coord[0] = zeta[0];
          for (unsigned j = 0; j < n_dim; j++)
          {
            node_coord[j + 1] = inode_pt->x(j);
          }
          set_node_coord.insert(node_coord);
        }


        if (travel_forward == true)
        {
          for (std::set<Vector<double>>::iterator it = set_node_coord.begin();
               it != set_node_coord.end();
               it++)
          {
            // Get the node coordinates
            Vector<double> node_coord = (*it);

            // Output the node coordinates
            const unsigned n_dim = node_coord.size() - 1;
            for (unsigned j = 0; j < n_dim; j++)
            {
              Boundary_output_file << node_coord[j + 1] << " ";
            }

            // Output the boundary coordinate
            Boundary_output_file << node_coord[0] << std::endl;

            Problem_Parameter::Ordered_bubbles[bubble][bubble_node].resize(2);
            Problem_Parameter::Ordered_bubbles[bubble][bubble_node][0] =
              node_coord[1];
            Problem_Parameter::Ordered_bubbles[bubble][bubble_node][1] =
              node_coord[2];

            bubble_node++;
          }
        }
        else
        {
          for (std::set<Vector<double>>::reverse_iterator it =
                 set_node_coord.rbegin();
               it != set_node_coord.rend();
               it++)
          {
            // Get the node coordinates
            Vector<double> node_coord = (*it);

            // Output the node coordinates
            const unsigned n_dim = node_coord.size() - 1;
            for (unsigned j = 0; j < n_dim; j++)
            {
              Boundary_output_file << node_coord[j + 1] << " ";
            }

            // Output the boundary coordinate
            Boundary_output_file << node_coord[0] << std::endl;

            Problem_Parameter::Ordered_bubbles[bubble][bubble_node].resize(2);
            Problem_Parameter::Ordered_bubbles[bubble][bubble_node][0] =
              node_coord[1];
            Problem_Parameter::Ordered_bubbles[bubble][bubble_node][1] =
              node_coord[2];

            bubble_node++;
          }
        }
      }
    }

    Boundary_output_file.close();


    /// Check the orientation of all bubbles
    for (unsigned bubble = 0; bubble < Problem_Parameter::N_Bubble; bubble++)
    {
      double area =
        get_polygon_area(Problem_Parameter::Ordered_bubbles[bubble]);

      // std::cout << "Bubble " << bubble  << " has area " << area << std::endl;

      if (area < 0)
      {
        std::reverse(Problem_Parameter::Ordered_bubbles[bubble].begin(),
                     Problem_Parameter::Ordered_bubbles[bubble].end());
      }
    }

    // std::cout << "Check neighbouring points " << std::endl;
    /// Remove any points which are too close to their neighbours
    /// This should remove only repeated nodes
    double distance_tol = 1e-8;
    for (unsigned bubble = 0; bubble < N_Bubble_initial; bubble++)
    {
      unsigned n_node = Problem_Parameter::Ordered_bubbles[bubble].size();
      double x0 = Problem_Parameter::Ordered_bubbles[bubble][0][0];
      double y0 = Problem_Parameter::Ordered_bubbles[bubble][0][1];
      Vector<Vector<double>> Accepted_vertices;
      Accepted_vertices.push_back(
        Problem_Parameter::Ordered_bubbles[bubble][0]);
      for (unsigned i_vertex = 1; i_vertex < n_node; i_vertex++)
      {
        double x1 = Problem_Parameter::Ordered_bubbles[bubble][i_vertex][0];
        double y1 = Problem_Parameter::Ordered_bubbles[bubble][i_vertex][1];
        if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) >
            distance_tol * distance_tol)
        {
          Accepted_vertices.push_back(
            Problem_Parameter::Ordered_bubbles[bubble][i_vertex]);
        }
        else
        {
          // std::cout << "Reject vertex at " <<
          // Problem_Parameter::Ordered_bubbles[bubble][i_vertex][0] << " " <<
          // Problem_Parameter::Ordered_bubbles[bubble][i_vertex][1] <<
          // std::endl;
        }
        x0 = x1;
        y0 = y1;
      }

      /// Check that first and last vertices are the same
      double x1 = Problem_Parameter::Ordered_bubbles[bubble][0][0];
      double y1 = Problem_Parameter::Ordered_bubbles[bubble][0][1];
      if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) >
          distance_tol * distance_tol)
      {
        Accepted_vertices.push_back(
          Problem_Parameter::Ordered_bubbles[bubble][0]);
        // std::cout << "Push back to close loop" << std::endl;
      }

      // std::cout << "Bubble " << bubble << " had size " <<
      // Problem_Parameter::Ordered_bubbles[bubble].size() << " and now has size
      // " << Accepted_vertices.size() << std::endl;
      Problem_Parameter::Ordered_bubbles[bubble] = Accepted_vertices;
    }
  }


  Data* Q_inv_data_pt;
  Data* St_data_pt;
  Data* Alpha_data_pt;
  Data* Drop_viscosity_data_pt;
  Data* Obstacle_height_data_pt;
  Data* Obstacle_width_data_pt;
  Data* G_data_pt;
  Data* U_data_pt;
  Data* Asymmetry_data_pt;
  Data* CoM_Y_data_pt;

  unsigned Number_of_negative_eigenvalues;

private:
  void create_free_surface_elements();

  void delete_free_surface_elements()
  {
    unsigned n_element = Free_surface_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      delete Free_surface_mesh_pt->element_pt(e);
    }
    Free_surface_mesh_pt->flush_element_and_node_storage();
  }

  void create_volume_constraint_elements();
  void delete_volume_constraint_elements()
  {
    std::cout << "Delete vector of constraint elements " << std::endl;
    for (unsigned i = 0; i < Vol_constraint_el_vector_pt.size(); i++)
    {
      delete Vol_constraint_el_vector_pt[i];
    }
    Vol_constraint_el_vector_pt.resize(0);
    //	for (unsigned i_bubble=0;i_bubble<Problem_Parameter::N_Bubble;
    //i_bubble++)
    //	{
    //		delete Vol_constraint_el_vector_pt[i_bubble];
    //	}

    std::cout << "Delete the individual elements " << std::endl;
    unsigned n_element = Volume_constraint_mesh_pt->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      delete Volume_constraint_mesh_pt->element_pt(e);
    }
    std::cout << "Flush the mesh " << std::endl;
    Volume_constraint_mesh_pt->flush_element_and_node_storage();
    std::cout << "Flush the equation mesh " << std::endl;

    Volume_constraint_mesh_equation_pt->flush_element_and_node_storage();

    std::cout << "Done deleting .. " << std::endl;
  }

  Mesh* CoM_X_constraint_mesh_pt;
  Mesh* CoM_Y_constraint_mesh_pt;
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

  /// Pointer to mesh of free surface elements
  Mesh* Free_surface_mesh_pt;

  /// Pointer to Fluid_mesh
  RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;

  Mesh* Volume_constraint_mesh_pt;
  Mesh* Volume_constraint_mesh_equation_pt;

  /// Vector storing pointer to the bubble polygons
  Vector<TriangleMeshPolygon*> Bubble_polygon_pt;

  /// Triangle mesh polygon for outer boundary
  TriangleMeshPolygon* Outer_boundary_polyline_pt;

  Vector<TriangleMeshPolygon*> Obstacle_polygon_pt;

  /// Pointer to a global bubble pressure datum


  Mesh* Inflow_mesh_pt;

  Vector<VolumeConstraintElement*> Vol_constraint_el_vector_pt;
  VolumeConstraintElement* CoM_X_constraint_el_pt;
  SelfReferentialVolumeConstraintElement* CoM_Y_constraint_el_pt;

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


}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
BubbleInChannelProblem<ELEMENT>::BubbleInChannelProblem()
{
  Constraints_are_steady = false;

  // std::cout << "Restarting constructor" << std::endl;

  // Allocate the timestepper -- this constructs the Problem's
  // time object with a sufficient amount of storage to store the
  // previous timsteps.
  bool adaptive_timestepping = true;
  add_time_stepper_pt(new BDF<2>(adaptive_timestepping));


  double initial_height = 0.00;
  double initial_width = 0.25;
  double initial_alpha = 10;
  double initial_drop_viscosity = 0.01;
  double initial_St = 1;
  double initial_G = 0;
  double initial_Q_inv = 20;
  double initial_Frame_speed = 0;
  double initial_asymmetry = 0.00;
  double initial_CoM_Y = 0.0;

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

  Drop_viscosity_data_pt = new Data(1);
  this->add_global_data(Drop_viscosity_data_pt);
  Drop_viscosity_data_pt->set_value(0, initial_drop_viscosity);
  Problem_Parameter::global_drop_viscosity_pt =
    Drop_viscosity_data_pt->value_pt(0);
  Drop_viscosity_data_pt->pin(0);


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


  U_data_pt = new Data(1);
  this->add_global_data(U_data_pt);
  U_data_pt->set_value(0, initial_Frame_speed);
  Problem_Parameter::global_Frame_speed_pt = U_data_pt->value_pt(0);
  U_data_pt->pin(0);


  unsigned index_of_traded_data = 0;
  CoM_X_constraint_el_pt = new VolumeConstraintElement(
    &Problem_Parameter::Centre_of_mass, U_data_pt, index_of_traded_data);

  CoM_Y_constraint_el_pt = new SelfReferentialVolumeConstraintElement(
    Problem_Parameter::global_CoM_Y_pt, CoM_Y_data_pt, index_of_traded_data);


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
  vertex_coord[0][0] = Problem_Parameter::Channel_start;
  vertex_coord[0][1] = -1.0;
  vertex_coord[1][0] = Problem_Parameter::Channel_start;
  vertex_coord[1][1] = 1.0;

  // Build the 1st boundary polyline
  boundary_polyline_pt[0] =
    new TriangleMeshPolyLine(vertex_coord, Inflow_boundary_id);

  // Second boundary polyline: Upper wall
  vertex_coord[0][0] = Problem_Parameter::Channel_start;
  vertex_coord[0][1] = 1.0;
  vertex_coord[1][0] = Problem_Parameter::Channel_end;
  vertex_coord[1][1] = 1.0;

  // Build the 2nd boundary polyline
  boundary_polyline_pt[1] =
    new TriangleMeshPolyLine(vertex_coord, Upper_wall_boundary_id);

  // Third boundary polyline: Outflow
  vertex_coord[0][0] = Problem_Parameter::Channel_end;
  vertex_coord[0][1] = 1.0;
  vertex_coord[1][0] = Problem_Parameter::Channel_end;
  vertex_coord[1][1] = -1.0;

  // Build the 3rd boundary polyline
  boundary_polyline_pt[2] =
    new TriangleMeshPolyLine(vertex_coord, Outflow_boundary_id);

  // Fourth boundary polyline: Bottom wall
  vertex_coord[0][0] = Problem_Parameter::Channel_end;
  vertex_coord[0][1] = -1.0;
  vertex_coord[1][0] = Problem_Parameter::Channel_start;
  vertex_coord[1][1] = -1.0;

  // Build the 4th boundary polyline
  boundary_polyline_pt[3] =
    new TriangleMeshPolyLine(vertex_coord, Bottom_wall_boundary_id);

  // Create the triangle mesh polygon for outer boundary
  Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);


  // Now define initial shape of bubble(s) with polygon
  //---------------------------------------------------

  if (Problem_Parameter::Reload_from_vector == true)
  {
    Problem_Parameter::Read_from_file = false;
  }

  Bubble_polygon_pt.resize(Problem_Parameter::N_Bubble);

  /// This parameter feeds into the interface remeshing. However, it's not clear
  /// (I can't remember) whether it refers to a physical or scaled length.
  double maximum_length = 0.02;
  Vector<Vector<double>> Bubble_center_vector;

  Bubble_center_vector.resize(Problem_Parameter::N_Bubble);
  if (Problem_Parameter::Reload_from_vector == true)
  {
    //	std::cout << "Reloading from vector " << std::endl;
    for (unsigned i_circle = 0; i_circle < Problem_Parameter::N_Bubble;
         i_circle++)
    {
      unsigned npoints_read =
        Problem_Parameter::All_the_bubbles[i_circle].size();
      //		std::cout << "Bubble " << i_circle << " has " << npoints_read << "
      //points" << std::endl;

      Vector<double> x_position(npoints_read, 0.0);
      Vector<double> y_position(npoints_read, 0.0);

      for (unsigned i = 0; i < npoints_read; i++)
      {
        x_position[i] = Problem_Parameter::All_the_bubbles[i_circle][i][0];
        y_position[i] = Problem_Parameter::All_the_bubbles[i_circle][i][1];
        //	std::cout << x_position[i] << " " << y_position[i] << std::endl;
      }
      Vector<double> Point_inside(2, 0.0);
      UnstructuredTwoDMeshGeometryBase::find_point_inside_polygon_helper(
        Problem_Parameter::All_the_bubbles[i_circle], Point_inside);
      double x_center = Point_inside[0];
      double y_center = Point_inside[1];

      Bubble_center_vector[i_circle] = Point_inside;

      Vector<TriangleMeshCurveSection*> bubble_polyline_pt(2);

      unsigned npoints = npoints_read / 2;

      Vector<Vector<double>> bubble_vertex(npoints);
      unsigned n = 0;

      //	std::cout << "Place " << npoints << " points into first boundary " <<
      //std::endl;
      for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
      {
        bubble_vertex[ipoint].resize(2);

        bubble_vertex[ipoint][0] = x_position[n];
        bubble_vertex[ipoint][1] = y_position[n];

        //			std::cout << bubble_vertex[ipoint][0] << " " <<
        //bubble_vertex[ipoint][1] << std::endl;
        n++;
      }
      bubble_polyline_pt[0] = new TriangleMeshPolyLine(
        bubble_vertex, First_bubble_boundary_id + 2 * i_circle);
      double first_x = bubble_vertex[0][0];
      double first_y = bubble_vertex[0][1];

      n = n - 1;

      npoints = npoints_read - npoints + 1;

      //		std::cout << "Place " << npoints << " points into second boundary "
      //<< std::endl;
      bubble_vertex.resize(npoints);
      for (unsigned ipoint = 0; ipoint < npoints; ipoint++)
      {
        bubble_vertex[ipoint].resize(2);

        //                std::cout << "Using entry " << n << " of up to " <<
        //                npoints_read-1 << std::endl;
        bubble_vertex[ipoint][0] = x_position[n];
        bubble_vertex[ipoint][1] = y_position[n];
        //			std::cout << bubble_vertex[ipoint][0] << " " <<
        //bubble_vertex[ipoint][1] << std::endl;

        n++;
      }


      /// Check loop is closed...
      double last_x = bubble_vertex[npoints - 1][0];
      double last_y = bubble_vertex[npoints - 1][1];
      if ((first_x - last_x) * (first_x - last_x) +
            (first_y - last_y) * (first_y - last_y) >
          1e-12)
      {
        std::cout << "There may be a mismatch between first and last point on "
                     "this bubble"
                  << std::endl;
        std::cout << "Will manually close the loop (in constructor) "
                  << std::endl;
        Vector<double> additional_point(2, 0.0);
        additional_point[0] = first_x;
        additional_point[1] = first_y;
        bubble_vertex.push_back(additional_point);
      }

      bubble_polyline_pt[1] = new TriangleMeshPolyLine(
        bubble_vertex, Second_bubble_boundary_id + 2 * i_circle);

      // Define coordinates of a point inside the bubble
      Vector<double> bubble_center(2);
      bubble_center[0] = x_center;
      bubble_center[1] = y_center;

      /// From vector
      bubble_polyline_pt[0]->set_maximum_length(maximum_length);
      bubble_polyline_pt[1]->set_maximum_length(maximum_length);

      // Create closed polygon from two polylines
      Bubble_polygon_pt[i_circle] = new TriangleMeshPolygon(bubble_polyline_pt);

      /// These tolerances correspond to interface curvatures.

      Bubble_polygon_pt[i_circle]->set_polyline_refinement_tolerance(0.05);
      Bubble_polygon_pt[i_circle]->set_polyline_unrefinement_tolerance(0.02);
      Bubble_polygon_pt[i_circle]
        ->enable_redistribution_of_segments_between_polylines();
      //	std::cout << "Bubble " << i_circle << std::endl;
    }
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
  double uniform_element_area = 0.4;

  // Use the TriangleMeshParameters object for gathering all
  // the necessary arguments for the TriangleMesh object
  TriangleMeshParameters triangle_mesh_parameters(outer_closed_curve_pt);

  // Define the holes on the boundary
  triangle_mesh_parameters.internal_closed_curve_pt() = bubble_closed_curve_pt;

  // Region coordinates is done in a loop over each drop/bubble.
  for (unsigned i_bubble = 0; i_bubble < Problem_Parameter::N_Bubble;
       i_bubble++)
  {
    //	std::cout << "Setting region " << i_bubble+1 << std::endl;
    triangle_mesh_parameters.add_region_coordinates(
      i_bubble + 1, Bubble_center_vector[i_bubble]);
    //	std::cout << "Drop " << i_bubble << " has center: " <<
    //Bubble_center_vector[i_bubble][0] << " " <<
    //Bubble_center_vector[i_bubble][1] << std::endl;
  }

  // Define the maximum element areas
  triangle_mesh_parameters.element_area() = uniform_element_area;

  // Create the mesh
  Fluid_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
    triangle_mesh_parameters, this->time_stepper_pt());

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Fluid_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  Fluid_mesh_pt->max_permitted_error() = 2e-3;
  Fluid_mesh_pt->min_permitted_error() = 5e-6;
  Fluid_mesh_pt->max_element_size() = 0.2;
  Fluid_mesh_pt->min_element_size() = 5e-3;
  //	Fluid_mesh_pt->disable_projection();
  //	Fluid_mesh_pt->disable_iterative_solver_for_projection();


  complete_problem_setup();

  Free_surface_mesh_pt = new Mesh;
  create_free_surface_elements();

  Problem_Parameter::Bubble_volume_vector.resize(Problem_Parameter::N_Bubble);
  Volume_constraint_mesh_pt = new Mesh;
  Volume_constraint_mesh_equation_pt = new Mesh;
  // create_volume_constraint_elements();
  // delete_volume_constraint_elements();
  // create_volume_constraint_elements();

  CoM_X_constraint_mesh_pt = new Mesh;
  create_CoM_X_constraint_elements();

  CoM_Y_constraint_mesh_pt = new Mesh;
  create_CoM_Y_constraint_elements();

  Inflow_mesh_pt = new Mesh;
  create_inflow_elements();

  this->add_sub_mesh(Fluid_mesh_pt);

  this->add_sub_mesh(this->Free_surface_mesh_pt);
  // this->add_sub_mesh(this->Volume_constraint_mesh_pt);
  // this->add_sub_mesh(this->Volume_constraint_mesh_equation_pt);
  this->add_sub_mesh(this->CoM_X_constraint_mesh_pt);
  this->add_sub_mesh(this->CoM_Y_constraint_mesh_pt);
  this->add_sub_mesh(this->Inflow_mesh_pt);

  this->build_global_mesh();

  // cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
  redo_equation_numbering();


  //	Use_finite_differences_for_continuation_derivatives = true;

  // linear_solver_pt()=new FD_LU;

#ifdef OOMPH_HAS_TRILINOS
  eigen_solver_pt() = new ANASAZI;
  static_cast<ANASAZI*>(eigen_solver_pt())->set_shift(-10);
#endif

  Max_residuals = 1e6;
  Max_newton_iterations = 20;
} // end_of_constructor


//============start_of_create_free_surface_elements======================
/// Create elements that impose the kinematic and dynamic bcs
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_free_surface_elements()
{
  //	std::cout << "Create free surface elements, in region 0 (bulk) only" <<
  //std::endl;
  unsigned region = 0;
  for (unsigned bubble = 0; bubble < Problem_Parameter::N_Bubble; bubble++)
  {
    for (unsigned b = First_bubble_boundary_id + 2 * bubble;
         b < Second_bubble_boundary_id + 1 + 2 * bubble;
         b++)
    {
      //		std::cout << "Setting free surface elements on boundary " << b <<
      //std::endl;
      unsigned n_element =
        Fluid_mesh_pt->nboundary_element_in_region(b, region);
      //		std::cout << "This is : " << n_element << " bulk elements"
      //<<std::endl;
      for (unsigned e = 0; e < n_element; e++)
      {
        ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
          Fluid_mesh_pt->boundary_element_in_region_pt(b, region, e));

        int face_index =
          Fluid_mesh_pt->face_index_at_boundary_in_region(b, region, e);

        HeleShawInterfaceElement<ELEMENT>* el_pt =
          new HeleShawInterfaceElement<ELEMENT>(bulk_elem_pt, face_index);

        Free_surface_mesh_pt->add_element_pt(el_pt);

        el_pt->set_boundary_number_in_bulk_mesh(b);

        el_pt->q_inv_pt() = Q_inv_data_pt->value_pt(0);

        el_pt->st_pt() = St_data_pt->value_pt(0);

        el_pt->alpha_pt() = Alpha_data_pt->value_pt(0);
        el_pt->drop_viscosity_pt() = Drop_viscosity_data_pt->value_pt(0);

        el_pt->upper_wall_fct_pt() = Problem_Parameter::upper_wall_fct;
        el_pt->drop_pressure_fct_pt() =
          Problem_Parameter::drop_pressure_function;

        el_pt->wall_speed_fct_pt() = Problem_Parameter::frame_speed_fct;

        el_pt->add_external_data(Q_inv_data_pt, true);
        el_pt->add_external_data(U_data_pt, true);

        el_pt->add_external_data(Alpha_data_pt, true);
        el_pt->add_external_data(Drop_viscosity_data_pt, true);
      }
    }
  }

  //	std::cout << "Finished create free surface " << std::endl;
}
// end of create_free_surface_elements

//============start_of_create_volume_constraint_elements=================
/// Create elements that impose volume constraint on the bubble
//=======================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_volume_constraint_elements()
{
  //	std::cout << "Create volume constraint " << std::endl;
  Vol_constraint_el_vector_pt.resize(Problem_Parameter::N_Bubble);
  for (unsigned i_bubble = 0; i_bubble < Problem_Parameter::N_Bubble;
       i_bubble++)
  {
    unsigned constraint_region = 1;
    /// Attach constraint to first node in first element in region 1.
    /// This should trade the value of the pressure there for the requirement
    /// that the volume matches.
    unsigned constraint_elt = 0;
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      Fluid_mesh_pt->region_element_pt(constraint_region, constraint_elt));
    Vol_constraint_el_vector_pt[i_bubble] = new VolumeConstraintElement(
      &Problem_Parameter::Bubble_volume_vector[i_bubble],
      el_pt->node_pt(0),
      constraint_region);
    Volume_constraint_mesh_equation_pt->add_element_pt(
      Vol_constraint_el_vector_pt[i_bubble]);
  }


  unsigned region = 0;
  for (unsigned bubble = 0; bubble < Problem_Parameter::N_Bubble; bubble++)
  {
    for (unsigned b = First_bubble_boundary_id + 2 * bubble;
         b < Second_bubble_boundary_id + 1 + 2 * bubble;
         b++)
    {
      unsigned n_element =
        Fluid_mesh_pt->nboundary_element_in_region(b, region);
      for (unsigned e = 0; e < n_element; e++)
      {
        ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
          Fluid_mesh_pt->boundary_element_in_region_pt(b, region, e));

        int face_index =
          Fluid_mesh_pt->face_index_at_boundary_in_region(b, region, e);

        HeleShawVolumeConstraintElement<ELEMENT>* el_pt =
          new HeleShawVolumeConstraintElement<ELEMENT>(bulk_elem_pt,
                                                       face_index);

        el_pt->set_volume_constraint_element(
          Vol_constraint_el_vector_pt[bubble]);
        el_pt->i_bubble = bubble;

        Volume_constraint_mesh_pt->add_element_pt(el_pt);
      }
    }
  }
  Constraints_are_steady = true;
  //	std::cout << "Finished create volume constraint elements " << std::endl;
}
// end of create_volume_constraint_elements


template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_CoM_X_constraint_elements()
{
  unsigned region = 0;
  CoM_X_constraint_mesh_pt->add_element_pt(CoM_X_constraint_el_pt);

  for (unsigned bubble = 0; bubble < Problem_Parameter::N_Bubble; bubble++)
  {
    for (unsigned b = First_bubble_boundary_id + 2 * bubble;
         b < Second_bubble_boundary_id + 1 + 2 * bubble;
         b++)
    {
      //		std::cout << "Setting CoM X constraint elements on boundary " << b
      //<< std::endl;
      unsigned n_element =
        Fluid_mesh_pt->nboundary_element_in_region(b, region);

      for (unsigned e = 0; e < n_element; e++)
      {
        ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
          Fluid_mesh_pt->boundary_element_in_region_pt(b, region, e));

        int face_index =
          Fluid_mesh_pt->face_index_at_boundary_in_region(b, region, e);

        HeleShawCoMXConstraintElement<ELEMENT>* el_pt =
          new HeleShawCoMXConstraintElement<ELEMENT>(bulk_elem_pt, face_index);

        el_pt->set_volume_constraint_element(CoM_X_constraint_el_pt);
        CoM_X_constraint_mesh_pt->add_element_pt(el_pt);
      }
    }
  }
  //	std::cout << "Finished create CoM X constraint elements " << std::endl;
}
// end of create_CoM_X_constraint_elements

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_CoM_Y_constraint_elements()
{
  unsigned region = 0;
  // Add volume constraint element to the mesh
  CoM_Y_constraint_mesh_pt->add_element_pt(CoM_Y_constraint_el_pt);


  for (unsigned b = First_bubble_boundary_id; b < Second_bubble_boundary_id + 1;
       b++)
  {
    //		std::cout << "Setting CoM Y constraint elements on boundary " << b <<
    //std::endl;
    // How many bulk fluid elements are adjacent to boundary b?
    unsigned n_element = Fluid_mesh_pt->nboundary_element_in_region(b, region);

    // Loop over the bulk fluid elements adjacent to boundary b?
    for (unsigned e = 0; e < n_element; e++)
    {
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_in_region_pt(b, region, e));

      int face_index =
        Fluid_mesh_pt->face_index_at_boundary_in_region(b, region, e);

      // Create new element
      HeleShawCoMYConstraintElement<ELEMENT>* el_pt =
        new HeleShawCoMYConstraintElement<ELEMENT>(bulk_elem_pt, face_index);

      // Set the "master" volume constraint element
      el_pt->set_volume_constraint_element(CoM_Y_constraint_el_pt);

      // Add it to the mesh
      CoM_Y_constraint_mesh_pt->add_element_pt(el_pt);
    }
  }

  //	std::cout << "Finished create CoM Y constraint elements " << std::endl;
}

//==start_of_complete_problem_setup=======================================
/// Set boundary conditions and complete the build of all elements
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::complete_problem_setup()
{
  // std::cout << "Complete problem setup:" << std::endl;

  map<unsigned, bool> is_on_bubble_bound;
  unsigned nbubble = Problem_Parameter::N_Bubble;
  for (unsigned ibubble = 0; ibubble < nbubble; ibubble++)
  {
    Vector<unsigned> bubble_bound_id =
      this->Bubble_polygon_pt[ibubble]->polygon_boundary_id();
    unsigned nbound = bubble_bound_id.size();
    for (unsigned ibound = 0; ibound < nbound; ibound++)
    {
      unsigned bound_id = bubble_bound_id[ibound];
      is_on_bubble_bound[bound_id] = true;
    }
  }

  // std::cout << "Loop over all nodes in fluid mesh, pinning both pressure
  // nodes" <<std::endl;
  const unsigned n_node = Fluid_mesh_pt->nnode();
  for (unsigned n = 0; n < n_node; n++)
  {
    Fluid_mesh_pt->node_pt(n)->set_value(0, 0.0);
    Fluid_mesh_pt->node_pt(n)->pin(0);
    Fluid_mesh_pt->node_pt(n)->set_value(1, 0.0);
    Fluid_mesh_pt->node_pt(n)->pin(1);
  }

  // std::cout << "Loop over all bulk elements, setting functions" <<std::endl;
  unsigned n_element = Fluid_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
    el_pt->constitutive_law_pt() = Problem_Parameter::Constitutive_law_pt;
    el_pt->upper_wall_fct_pt() = Problem_Parameter::upper_wall_fct;
  }

  // std::cout << "Loop over all nodes in region 0, unpinning pressure 0, which
  // is the relevant pressure" <<std::endl;
  unsigned n_region_element = Fluid_mesh_pt->nregion_element(0);
  for (unsigned e = 0; e < n_region_element; e++)
  {
    ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(Fluid_mesh_pt->region_element_pt(0, e));
    el_pt->set_pressure_node_value(0);
    unsigned nodes_in_element = el_pt->nnode();
    for (unsigned i_node = 0; i_node < nodes_in_element; i_node++)
    {
      el_pt->node_pt(i_node)->unpin(0);
    }
    el_pt->Viscosity = 0;
  }

  // std::cout << "Loop over all nodes in region 1, unpinning pressure 1, which
  // is the relevant pressure" <<std::endl;
  for (unsigned i_bubble = 0; i_bubble < Problem_Parameter::N_Bubble;
       i_bubble++)
  {
    unsigned drop_region = i_bubble + 1;
    n_region_element = Fluid_mesh_pt->nregion_element(drop_region);
    for (unsigned e = 0; e < n_region_element; e++)
    {
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->region_element_pt(drop_region, e));
      el_pt->set_pressure_node_value(1);
      unsigned nodes_in_element = el_pt->nnode();
      for (unsigned i_node = 0; i_node < nodes_in_element; i_node++)
      {
        el_pt->node_pt(i_node)->unpin(1);
      }
      el_pt->Viscosity = 1;
    }
  }

  /// Loop over nodes on solid boundary to pin position, and pressure at
  /// outflow.
  unsigned nbound = Fluid_mesh_pt->nboundary();
  for (unsigned ibound = 0; ibound < nbound; ibound++)
  {
    unsigned num_nod = Fluid_mesh_pt->nboundary_node(ibound);
    for (unsigned inod = 0; inod < num_nod; inod++)
    {
      Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound, inod);
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
      }
    }
  }

} // end of complete_problem_setup

//========================================================================
/// Compute error estimates and assign to elements for plotting
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::compute_error_estimate(double& max_err,
                                                             double& min_err)
{
  //	std::cout << "Compute error estimate" << std::endl;
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

  //	std::cout << "Max error is " << max_err << std::endl;
  //	std::cout << "Min error is " << min_err << std::endl;
}

//========================================================================
/// Compute error estimates and assign to elements for plotting
//========================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::solve_for_eigenproblem()
{
  bool doc_eigenvector =
    false; /// Will doc perturbed solutions for 10 eigenvectors.
  bool retain_perturbed_solution_at_function_exit =
    false; /// System state at function
  /// exit will be perturbed by eigenvector 0. Good for finite amplitude
  /// calculations.

  /// If we want to output solutions, we should backup the initial state.
  unsigned n_dof = ndof();
  Vector<double> backup_1(n_dof);
  // Store original state.
  for (unsigned n = 0; n < n_dof; n++)
  {
    backup_1[n] = dof(n);
  }

  if (doc_eigenvector == true)
  {
    //		std::cout << "Doc initial solution" << std::endl;
    doc_solution();
  }
  Problem_Parameter::vector_of_eigenvalues_rp.resize(10);
  Problem_Parameter::vector_of_eigenvalues_ip.resize(10);

  /// Reset eigenvalue guesses.
  for (unsigned n = 0; n < 10; n++)
  {
    Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
    Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
  }

  // Set external storage for the eigenvalues
  Vector<complex<double>> eigenvalues;
  Vector<DoubleVector> eigenvectors;

  // Desired number of eigenvalues
  unsigned n_eval = 10;
  int sensible_eigenvalues = 0;
  // Solve the eigenproblem
  std::cout << "Now attempt to solve eigenproblem, n_eval = " << n_eval
            << std::endl;
  solve_eigenproblem(n_eval, eigenvalues, eigenvectors);
  std::cout << "N_eval is " << n_eval << std::endl;
  std::cout << "Eigenvalues_size is " << eigenvalues.size() << std::endl;

  /// Describe eigenvalues
  std::cout << "And describe eigenvalues" << std::endl;
  n_eval = eigenvalues.size();
  Number_of_negative_eigenvalues = 0;
  for (unsigned n = 0; n < n_eval; n++)
  {
    if (isinf(real(eigenvalues[n])) != true &&
        isnan(real(eigenvalues[n])) != true)
    {
      std::cout << "Eigenvalue " << eigenvalues[n] << std::endl;
      sensible_eigenvalues++;
      if (real(eigenvalues[n]) < 0)
      {
        Number_of_negative_eigenvalues++;
      }
    }
    if (doc_eigenvector == true)
    {
      std::cout << "Doc eigenvector with eigenvalue: " << real(eigenvalues[n])
                << " " << imag(eigenvalues[n]) << " " << std::endl;
      /// The function assign_eigenvector_to_dofs(eigenvectors[n]) is not
      /// helpful for finding the perturbed interface position.

      /// Perturb original solution. Choice of factor 10 is arbitrary.
      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = backup_1[i] + 10 * eigenvectors[n][i];
      }
      actions_after_change_in_bifurcation_parameter();
      doc_solution();
    }
  }

  for (unsigned n = 0; n < n_dof; n++)
  {
    dof(n) = backup_1[n];
  }
  actions_after_change_in_bifurcation_parameter();

  for (unsigned n = 0; n < 10; n++)
  {
    if (n_eval >= n + 1)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = real(eigenvalues[n]);
      Problem_Parameter::vector_of_eigenvalues_ip[n] = imag(eigenvalues[n]);
    }
  }

  std::cout << "There are " << sensible_eigenvalues << " reasonable eigenvalues"
            << std::endl;

  if (retain_perturbed_solution_at_function_exit == true)
  {
    std::cout << "Keep perturbation with eigenvector 0" << std::endl;
    unsigned perturbation_evec = 0;

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) = backup_1[i] + 2 * eigenvectors[perturbation_evec][i];
    }
    actions_after_change_in_bifurcation_parameter();
    doc_solution();
  }
}

//============start_of_create_volume_constraint_elements=================
/// Create elements that impose volume constraint on the bubble
//=======================================================================
template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::assess_bubble_volumes(
  bool reset_bubble_volumes)
{
  Vector<double> Temp_volume_vector(Problem_Parameter::N_Bubble, 0.0);
  unsigned n_element = Volume_constraint_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    HeleShawVolumeConstraintElement<ELEMENT>* el_pt =
      dynamic_cast<HeleShawVolumeConstraintElement<ELEMENT>*>(
        Volume_constraint_mesh_pt->element_pt(e));
    unsigned i_b = el_pt->i_bubble;
    double dV = 0;
    el_pt->calculate_volume_contribution(dV);
    Temp_volume_vector[i_b] += dV;
  }

  for (unsigned i_b = 0; i_b < Problem_Parameter::N_Bubble; i_b++)
  {
    std::cout << "Total contribution to bubble " << i_b << " is "
              << Temp_volume_vector[i_b] << std::endl;
    std::cout << "Volume constraint " << i_b << " residual is "
              << Temp_volume_vector[i_b] -
                   Problem_Parameter::Bubble_volume_vector[i_b]
              << std::endl;
    if (reset_bubble_volumes == true)
    {
      Problem_Parameter::Bubble_volume_vector[i_b] = Temp_volume_vector[i_b];
    }
  }
}
