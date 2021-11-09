#ifndef HELE_SHAW_BUBBLE_PROBLEM_WITH_INTEGRALS_PERTURBATION_HEADER
#define HELE_SHAW_BUBBLE_PROBLEM_WITH_INTEGRALS_PERTURBATION_HEADER
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
    delete_CoM_X_constraint_for_Q_elements();
    delete_CoM_Y_constraint_elements();
    delete_inflow_elements();

    delete CoM_X_constraint_mesh_pt;
    delete CoM_X_constraint_for_Q_mesh_pt;
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
    delete_CoM_X_constraint_for_Q_elements();
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
    create_CoM_X_constraint_for_Q_elements();
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

  /// Update the problem specs before solve
  void actions_before_newton_solve() {}

  void jack_solve(double Q_unsteady,
                  bool jack_arc,
                  unsigned n_iter,
                  double sign);
  void dump_it(ofstream& dump_file);

  void restart(ifstream& restart_file);

  void set_initial_condition();
  void set_initial_condition_asymm();

  void test_jacobian_call();

  void read_eigenvalues_eigenvectors(complex<double>& evalue,
                                     Vector<complex<double>>& g,
                                     int mode)
  {
    int ndof = this->ndof();
    std::ifstream indata;
    char filename[100];

    vector<double> g_real(ndof);
    vector<double> g_imag(ndof);

    double evalue_real;
    double evalue_imag;

    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvectors_%i.dat").c_str(),
      mode);
    indata.open(filename);

    for (int i = 0; i < ndof; i++)
    {
      indata >> g_real[i];
    }

    indata.close();

    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvectors_%i.dat").c_str(),
      mode + 1);

    indata.open(filename);

    for (int i = 0; i < ndof; i++)
    {
      indata >> g_imag[i];
    }

    indata.close();

    for (int i = 0; i < ndof; i++)
    {
      g[i] = g_real[i] + 1i * g_imag[i];
    }

    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvalues_0.dat").c_str());
    indata.open(filename);

    indata >> evalue_real;
    indata >> evalue_imag;
    indata.close();

    evalue = evalue_real + 1i * evalue_imag;
  }

  void read_adjoint_eigenvalues_eigenvectors(complex<double>& evalue,
                                             Vector<complex<double>>& g,
                                             int mode)
  {
    int ndof = this->ndof();
    std::ifstream indata;
    char filename[100];

    vector<double> g_real(ndof);
    vector<double> g_imag(ndof);

    double evalue_real;
    double evalue_imag;

    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvectors_adjoint_%i.dat")
        .c_str(),
      mode);
    indata.open(filename);

    for (int i = 0; i < ndof; i++)
    {
      indata >> g_real[i];
    }

    indata.close();

    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvectors_adjoint_%i.dat")
        .c_str(),
      mode + 1);

    indata.open(filename);

    for (int i = 0; i < ndof; i++)
    {
      indata >> g_imag[i];
    }

    indata.close();

    for (int i = 0; i < ndof; i++)
    {
      g[i] = g_real[i] + 1i * g_imag[i];
    }

    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvalues_adjoint_0.dat")
        .c_str());
    indata.open(filename);

    indata >> evalue_real;
    indata >> evalue_imag;

    indata.close();

    evalue = evalue_real + 1i * evalue_imag;
  }

  void fill_in_steady_fourth_order_solution(Vector<complex<double>>& varphi_4,
                                            Vector<complex<double>>& varphi_3,
                                            int delta);

  void fill_in_second_order_solution(Vector<complex<double>>& g,
                                     double Q_c,
                                     double omega_c,
                                     int mode,
                                     Vector<complex<double>>& varphi_0,
                                     Vector<complex<double>>& varphi_2,
                                     Vector<complex<double>>& varphi_3);

  void fill_in_landau_coefficients(Vector<complex<double>>& g,
                                   Vector<complex<double>>& g_adjoint,
                                   Vector<complex<double>>& varphi_0,
                                   Vector<complex<double>>& varphi_2,
                                   Vector<complex<double>>& varphi_3,
                                   double Q_c,
                                   double omega_c);

  void mult_constant(CRDoubleMatrix& A, double mult);

  void mult_constant(CRComplexMatrix& A, complex<double> mult);

  void get_hessian(DoubleVector& g1,
                   DoubleVector& g2,
                   vector<double>& u0,
                   DoubleVector& f);

  void get_tressian(DoubleVector& g1,
                    DoubleVector& g2,
                    vector<double>& u0,
                    DoubleVector& f);

  void get_mass(DoubleVector& g1,
                DoubleVector& g2,
                vector<double>& u0,
                DoubleVector& f);

  void test_identity_matrix();

  void test_complex_identity_matrix();

  void test_complex_solve();

  void create_complex_mass_matrix(CRComplexMatrix& mass_complex,
                                  complex<double> mult);

  void create_complex_jacobian_matrix(CRComplexMatrix& mass_complex,
                                      complex<double> mult);

  void edge_track()
  {
    // Vary epsilon with theta held constant
    int N = this->ndof();
    // Perturbation of one side of orbit
    vector<double> v1(N);
    // Perturbation of other side of orbit
    vector<double> v2(N);
  }

  void assign_perturbation_time(Vector<double>& perturb)
  {
    int N = this->ndof();
    DoubleVector u(this->dof_distribution_pt());
    DoubleVector v(this->dof_distribution_pt());
    unsigned tt;

    for (int i = 0; i < N; i++)
    {
      u[i] = this->dof(i);
      v[i] = perturb[i];
    }

    tt = 0;
    set_dofs(tt, v);
    tt = 1;
    set_dofs(tt, u);
    tt = 2;
    set_dofs(tt, u);
  }

  void solve_for_eigenproblem_com()
  {
    bool doc_eigenvector = false;
    bool silent = false;

    Vector<complex<double>> eigenvalues;
    Vector<DoubleVector> eigenvectors;

    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    if (doc_eigenvector == true)
    {
      if (silent == false)
      {
        std::cout << "Doc Initial Solution" << std::endl;
      }

      doc_solution();
    }
    // Reset Eigenvalues

    for (unsigned n = 0; n < 10; n++)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
      Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
    }

    // Vector<complex<double> > eigenvalues;
    // Vector<DoubleVector> eigenvectors;

    unsigned n_eval = 10;
    int sensible_eigenvalues = 0;
    // Solve for eigenproblem
    if (silent == false)

    {
      std::cout << "Now attempt to solve the eigenproblem" << std::endl;
    }
    solve_eigenproblem(n_eval, eigenvalues, eigenvectors);

    if (silent == false)
    {
      std::cout << "N_eval is " << n_eval << std::endl;

      std::cout << "Eigenvalues size is" << eigenvalues.size() << std::endl;

      // Describe Eigenvalues
      std::cout << "And describe eigenvalues" << std::endl;
    }

    n_eval = eigenvalues.size();
    int Number_of_negative_eigenvalues = 0;
    for (unsigned n = 0; n < n_eval; n++)
    {
      if (isinf(real(eigenvalues[n])) != true &&
          isnan(real(eigenvalues[n])) != true)
      {
        if (silent == false)
        {
          std::cout << "Eigenvalue " << eigenvalues[n] << std::endl;
        }
        sensible_eigenvalues++;
        if (real(eigenvalues[n]) < 0)
        {
          Number_of_negative_eigenvalues++;
        }
      }
      if (doc_eigenvector == true)
      {
        if (silent == false)
        {
          std::cout << "Doc eigenvector with eigenvalue" << real(eigenvalues[n])
                    << " " << imag(eigenvalues[n]) << " " << std::endl;
        }

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

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    for (unsigned i = 0; i < n_dof; i++)
    {
      backup_1[i] = dof(i);
    }
    doc_solution();

    int eigenmode;
    int conjugate;
    double epsilon;

    std::cout << "=====================" << std::endl;
    std::cout << "What eigenmode? (0-9)" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cin >> eigenmode;
    std::cout << "=====================" << std::endl;
    std::cout << "What is the conjugate? (0-9)" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cin >> conjugate;
    std::cout << "What perturbation? (any double)" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cin >> epsilon;

    int N_theta = 12;
    int N_eps = 100;
    double dtheta = MathematicalConstants::Pi / 12.0;
    vector<double> theta(N_theta);
    vector<double> eps(N_eps);
    theta[0] = 0.0;
    eps[0] = 13.0;
    for (int i = 1; i < N_theta; i++)
    {
      theta[i] = theta[i - 1] + dtheta;
    }
    for (int i = 1; i < N_eps; i++)
    {
      eps[i] = eps[i - 1] + 1.0;
    }

    ofstream outdata;
    char filename[100];
    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenmode_data_com_y.dat")
        .c_str());
    outdata.open(filename);
    for (int l = 0; l < N_eps; l++)
    {
      for (int m = 0; m < N_theta; m++)
      {
        for (unsigned i = 0; i < n_dof; i++)
        {
          dof(i) =
            backup_1[i] + eps[l] * (cos(theta[m]) * eigenvectors[eigenmode][i] -
                                    sin(theta[m]) * eigenvectors[conjugate][i]);
        }
        std::cout << eps[l] << std::endl;
        std::cout << theta[m] << std::endl;
        output_eigenvector_perturbation(13 + l, m);
        //	      doc_solution();

        //      actions_after_change_in_bifurcation_parameter();
        outdata << eps[l] << " " << theta[m] << " " << get_CoM_Y() << " "
                << get_y_max() << " " << get_y_min() << std::endl;
      }
    }
    outdata.close();
    exit(1);

    for (int j = 0; j < N_theta; j++)
    {
      epsilon = 0.0;
      sprintf(filename,
              (Problem_Parameter::Doc_info.directory() +
               "eigenmode_data_eps_theta_%i.dat")
                .c_str(),
              j);
      outdata.open(filename);
      for (int m = 0; m < 30; m++)
      {
        for (unsigned i = 0; i < n_dof; i++)
        {
          dof(i) = backup_1[i] +
                   epsilon * (cos(theta[j]) * eigenvectors[eigenmode][i] -
                              sin(theta[j]) * eigenvectors[conjugate][i]);
        }
        actions_after_change_in_bifurcation_parameter();
        outdata << theta[j] << " " << get_CoM_Y() << " " << get_y_max() << " "
                << get_y_min() << " " << epsilon << std::endl;
        epsilon += 1.0;
      }
      outdata.close();
    }
  }

  double get_y_max()
  {
    int boundary = 4;
    unsigned nnodes = this->Fluid_mesh_pt->nboundary_node(boundary);
    Node* bound_node_pt = this->Fluid_mesh_pt->boundary_node_pt(boundary, 0);
    Vector<double> y(nnodes);
    double yy;
    //      double xx;

    for (unsigned i = 0; i < nnodes; i++)
    {
      bound_node_pt = this->Fluid_mesh_pt->boundary_node_pt(boundary, i);
      // xx = bound_node_pt->x(0); //x ordinate
      yy = bound_node_pt->x(1); // y ordinate
      y[i] = yy;
    }

    double y_max = y[0];
    for (unsigned i = 1; i < nnodes; i++)
    {
      if (y[i] > y_max)
      {
        y_max = y[i];
      }
      else
      {
      }
    }

    return y_max;
  }

  double get_y_min()
  {
    int boundary = 5;
    unsigned nnodes = this->Fluid_mesh_pt->nboundary_node(boundary);
    Node* bound_node_pt = this->Fluid_mesh_pt->boundary_node_pt(boundary, 0);
    Vector<double> y(nnodes);
    double yy;

    for (unsigned i = 0; i < nnodes; i++)
    {
      bound_node_pt = this->Fluid_mesh_pt->boundary_node_pt(boundary, i);
      yy = bound_node_pt->x(1);
      y[i] = yy;
    }

    double y_min = y[0];
    for (unsigned i = 1; i < nnodes; i++)
    {
      if (y[i] < y_min)
      {
        y_min = y[i];
      }
      else
      {
      }
    }

    return y_min;
  }


  void solve_for_eigenproblem()
  {
    Vector<complex<double>> eigenvalues;
    Vector<DoubleVector> eigenvectors;

    solve_for_eigenproblem(eigenvalues, eigenvectors);
  }

  void test_inner_product(vector<complex<double>>& g,
                          vector<complex<double>>& g_adj,
                          complex<double> evalue)
  {
    std::cout << "========================================" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Testing that the inner product is correct" << std::endl;

    int N = this->ndof();
    DoubleVector g_real(this->dof_distribution_pt());
    DoubleVector g_imag(this->dof_distribution_pt());
    DoubleVector g_adj_real(this->dof_distribution_pt());
    DoubleVector g_adj_imag(this->dof_distribution_pt());

    DoubleVector A(this->dof_distribution_pt());
    DoubleVector B(this->dof_distribution_pt());
    DoubleVector C(this->dof_distribution_pt());
    DoubleVector D(this->dof_distribution_pt());
    vector<complex<double>> result(N);

    CRDoubleMatrix J_real(this->dof_distribution_pt());
    CRDoubleMatrix J_imag(this->dof_distribution_pt());
    CRDoubleMatrix M_real(this->dof_distribution_pt());
    CRDoubleMatrix M_imag(this->dof_distribution_pt());

    for (int i = 0; i < N; i++)
    {
      g_real[i] = real(g[i]);
      g_imag[i] = imag(g[i]);
      g_adj_real[i] = real(g_adj[i]);
      g_adj_imag[i] = imag(g_adj[i]);
    }

    set_steady();

    Problem::get_eigenproblem_matrices(M_real, J_real);
    Problem::get_eigenproblem_matrices(M_imag, J_imag);

    M_real.multiply(g_real, A);
    M_imag.multiply(g_imag, B);

    J_real.multiply(g_real, C);
    J_imag.multiply(g_imag, D);

    double inn_prod1_real = 0.0;
    double inn_prod1_imag = 0.0;

    for (int i = 0; i < N; i++)
    {
      inn_prod1_real += A[i] * g_adj_real[i] - B[i] * g_adj_imag[i];
      inn_prod1_imag += B[i] * g_adj_real[i] + A[i] * g_adj_imag[i];
    }

    double inn_prod2_real = 0.0;
    double inn_prod2_imag = 0.0;

    for (int i = 0; i < N; i++)
    {
      inn_prod2_real += C[i] * g_adj_real[i] - D[i] * g_adj_imag[i];
      inn_prod2_imag += D[i] * g_adj_real[i] + C[i] * g_adj_imag[i];
    }

    complex<double> evalue1 = (inn_prod2_real + 1i * inn_prod2_imag) /
                              (inn_prod1_real + 1i * inn_prod1_imag);

    std::cout << " " << std::endl;
    std::cout << "Actual eigenvalue is " << evalue << std::endl;
    std::cout << "Inner product eigenvalue is " << evalue1 << std::endl;
    std::cout << " " << std::endl;
  }

  void test_eigenvector_eigenvalue(complex<double>& evalue,
                                   vector<complex<double>>& g,
                                   int mode)
  {
    std::cout << "========================================" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Testing that the eigenvector is correct" << std::endl;
    //      output_eigenmatrices();

    int m = Problem_Parameter::Doc_info.number();
    int N = this->ndof();
    DoubleVector g_real(this->dof_distribution_pt());
    DoubleVector g_imag(this->dof_distribution_pt());
    DoubleVector A(this->dof_distribution_pt());
    DoubleVector B(this->dof_distribution_pt());
    DoubleVector C(this->dof_distribution_pt());
    DoubleVector D(this->dof_distribution_pt());
    vector<complex<double>> result(N);

    double norm = 0.0;
    int index = 0;
    double max = 10000;

    CRDoubleMatrix J_real(this->dof_distribution_pt());
    CRDoubleMatrix J_imag(this->dof_distribution_pt());
    CRDoubleMatrix M_real(this->dof_distribution_pt());
    CRDoubleMatrix M_imag(this->dof_distribution_pt());
    double evalue_real = real(evalue);
    double evalue_imag = imag(evalue);

    for (int i = 0; i < N; i++)
    {
      g_real[i] = real(g[i]);
      g_imag[i] = imag(g[i]);
    }

    set_steady();

    Problem::get_eigenproblem_matrices(M_real, J_real);
    Problem::get_eigenproblem_matrices(M_imag, J_imag);

    J_real.multiply(g_real, A);
    J_imag.multiply(g_imag, B);

    M_real.multiply(g_real, C);
    M_imag.multiply(g_imag, D);

    for (int i = 0; i < N; i++)
    {
      result[i] = A[i] + 1i * B[i] -
                  (evalue_real + 1i * evalue_imag) * (C[i] + 1i * D[i]);
    }

    max = sqrt(pow(real(result[0]), 2) + pow(imag(result[0]), 2));

    for (int i = 0; i < N; i++)
    {
      if (sqrt(pow(real(result[i]), 2) + pow(imag(result[i]), 2)) > max)
      {
        index = i;
        max = sqrt(pow(real(result[i]), 2) + pow(imag(result[i]), 2));
      }

      norm += sqrt(pow(real(result[i]), 2) + pow(imag(result[i]), 2));
    }

    std::cout << "The norm of M*s*g - J*g is " << norm << std::endl;
    std::cout << " " << std::endl;
    std::cout << "The max of M*s*g - J*g is " << result[index] << " at "
              << index << std::endl;

    std::cout << " " << std::endl;

    ofstream outdata;
    char filename[100];
    sprintf(filename,
            (Problem_Parameter::Doc_info.directory() + "result_%i.dat").c_str(),
            m);

    outdata.open(filename);
    for (int i = 0; i < N; i++)
    {
      outdata << real(result[i]) << " " << imag(result[i]) << std::endl;
    }
    outdata.close();

    Problem_Parameter::Doc_info.number() += 1;
  }

  void test_adjoint_eigenvector_eigenvalue(complex<double>& evalue,
                                           vector<complex<double>>& g,
                                           int mode)
  {
    std::cout << "========================================" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Testing that the adjoint eigenvector is correct" << std::endl;

    //      output_eigenmatrices();

    int m = Problem_Parameter::Doc_info.number();
    int N = this->ndof();
    DoubleVector g_real(this->dof_distribution_pt());
    DoubleVector g_imag(this->dof_distribution_pt());
    DoubleVector A(this->dof_distribution_pt());
    DoubleVector B(this->dof_distribution_pt());
    DoubleVector C(this->dof_distribution_pt());
    DoubleVector D(this->dof_distribution_pt());
    vector<complex<double>> result(N);

    double norm = 0.0;
    int index = 0;
    double max = 10000;

    CRDoubleMatrix J_real(this->dof_distribution_pt());
    CRDoubleMatrix J_imag(this->dof_distribution_pt());
    CRDoubleMatrix M_real(this->dof_distribution_pt());
    CRDoubleMatrix M_imag(this->dof_distribution_pt());
    CRDoubleMatrix J_T_real(this->dof_distribution_pt());
    CRDoubleMatrix J_T_imag(this->dof_distribution_pt());
    CRDoubleMatrix M_T_real(this->dof_distribution_pt());
    CRDoubleMatrix M_T_imag(this->dof_distribution_pt());

    double evalue_real = real(evalue);
    double evalue_imag = imag(evalue);

    for (int i = 0; i < N; i++)
    {
      g_real[i] = real(g[i]);
      g_imag[i] = imag(g[i]);
    }

    set_steady();

    Problem::get_eigenproblem_matrices(M_real, J_real);
    Problem::get_eigenproblem_matrices(M_imag, J_imag);

    CRDoubleMatrix* J_T_real_pt = &J_T_real;
    CRDoubleMatrix* J_T_imag_pt = &J_T_imag;
    CRDoubleMatrix* M_T_real_pt = &M_T_real;
    CRDoubleMatrix* M_T_imag_pt = &M_T_imag;

    J_real.get_matrix_transpose(J_T_real_pt);
    J_imag.get_matrix_transpose(J_T_imag_pt);
    M_real.get_matrix_transpose(M_T_real_pt);
    M_imag.get_matrix_transpose(M_T_imag_pt);

    J_T_real.multiply(g_real, A);
    J_T_imag.multiply(g_imag, B);

    M_T_real.multiply(g_real, C);
    M_T_imag.multiply(g_imag, D);

    for (int i = 0; i < N; i++)
    {
      result[i] = A[i] + 1i * B[i] -
                  (evalue_real + 1i * evalue_imag) * (C[i] + 1i * D[i]);
    }

    max = sqrt(pow(real(result[0]), 2) + pow(imag(result[0]), 2));

    for (int i = 0; i < N; i++)
    {
      if (sqrt(pow(real(result[i]), 2) + pow(imag(result[i]), 2)) > max)
      {
        index = i;
        max = sqrt(pow(real(result[i]), 2) + pow(imag(result[i]), 2));
      }

      norm += sqrt(pow(real(result[i]), 2) + pow(imag(result[i]), 2));
    }

    std::cout << "The norm of M*s*g - J*g is " << norm << std::endl;
    std::cout << " " << std::endl;
    std::cout << "The max of M*s*g - J*g is " << result[index] << " at "
              << index << std::endl;
    std::cout << " " << std::endl;

    ofstream outdata;
    char filename[100];
    sprintf(filename,
            (Problem_Parameter::Doc_info.directory() + "result_%i.dat").c_str(),
            m);

    outdata.open(filename);
    for (int i = 0; i < N; i++)
    {
      outdata << real(result[i]) << " " << imag(result[i]) << std::endl;
    }
    outdata.close();

    Problem_Parameter::Doc_info.number() += 1;
  }

  void solve_for_eigenproblem_sum()
  {
    Vector<complex<double>> eigenvalues;
    Vector<DoubleVector> eigenvectors;

    solve_for_eigenproblem_sum(eigenvalues, eigenvectors);
  }


  void solve_for_eigenproblem_sum(Vector<complex<double>>& eigenvalues,
                                  Vector<DoubleVector>& eigenvectors)

  {
    bool doc_eigenvector = false;
    bool silent = false;

    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    if (doc_eigenvector == true)
    {
      if (silent == false)
      {
        std::cout << "Doc Initial Solution" << std::endl;
      }

      doc_solution();
    }
    // Reset Eigenvalues

    for (unsigned n = 0; n < 10; n++)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
      Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
    }

    // Vector<complex<double> > eigenvalues;
    // Vector<DoubleVector> eigenvectors;

    unsigned n_eval = 10;
    int sensible_eigenvalues = 0;
    // Solve for eigenproblem
    if (silent == false)
    {
      std::cout << "Now attempt to solve the eigenproblem" << std::endl;
    }

    solve_eigenproblem(n_eval, eigenvalues, eigenvectors);

    if (silent == false)
    {
      std::cout << "N_eval is " << n_eval << std::endl;

      std::cout << "Eigenvalues size is" << eigenvalues.size() << std::endl;

      // Describe Eigenvalues
      std::cout << "And describe eigenvalues" << std::endl;
    }

    n_eval = eigenvalues.size();
    int Number_of_negative_eigenvalues = 0;
    for (unsigned n = 0; n < n_eval; n++)
    {
      if (isinf(real(eigenvalues[n])) != true &&
          isnan(real(eigenvalues[n])) != true)
      {
        if (silent == false)
        {
          std::cout << "Eigenvalue " << eigenvalues[n] << std::endl;
        }
        sensible_eigenvalues++;
        if (real(eigenvalues[n]) < 0)
        {
          Number_of_negative_eigenvalues++;
        }
      }
      if (doc_eigenvector == true)
      {
        if (silent == false)
        {
          std::cout << "Doc eigenvector with eigenvalue" << real(eigenvalues[n])
                    << " " << imag(eigenvalues[n]) << " " << std::endl;
        }

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

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    for (unsigned i = 0; i < n_dof; i++)
    {
      backup_1[i] = dof(i);
    }

    int eigenmode1;
    double epsilon1;
    double theta1;
    int conjugate1;
    int eigenmode2;
    double epsilon2;
    double theta2;
    int conjugate2;
    int answer = 0;

    for (unsigned m = 0; m < 1000; m++)
    {
      if (answer == 0)
      {
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "Choices, choices....." << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "What is your first eigenmode? (0-9)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> eigenmode1;
        std::cout << "What is the conjugate? (0-9)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> conjugate1;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "What phase? (as a multiple of pi/12)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> theta1;
        std::cout << "What perturbation? (any double)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> epsilon1;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "What is your second eigenmode? (0-9)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> eigenmode2;
        std::cout << "What is the conjugate? (0-9)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> conjugate2;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "What phase? (as a multiple of pi/12)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> theta2;
        std::cout << "What perturbation? (any double)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> epsilon2;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "Sure you're happy? (true = 1/false = 0)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> answer;
      }
      else
      {
        break;
      }
    }

    theta1 = (1.0 / 12.0) * theta1 * MathematicalConstants::Pi;
    theta2 = (1.0 / 12.0) * theta2 * MathematicalConstants::Pi;

    ofstream some_file;
    char filename[100];

    sprintf(filename, "perturbation.dat");

    some_file.open(filename);
    some_file << eigenmode1 << " " << epsilon1 << " "
              << " " << theta1 << endl;
    some_file << eigenmode2 << " " << epsilon2 << " "
              << " " << theta2 << endl;
    some_file.close();

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) =
        100 * (cos(theta1) * eigenvectors[eigenmode1][i] -
               sin(theta1) * eigenvectors[conjugate1][i]); // First eigenmode
      dof(i) +=
        100 * (cos(theta1) * eigenvectors[eigenmode2][i] -
               sin(theta2) * eigenvectors[conjugate2][i]); // First eigenmode
    }

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) =
        backup_1[i] + epsilon1 * (cos(theta1) * eigenvectors[eigenmode1][i] -
                                  sin(theta1) * eigenvectors[conjugate1][i]);
      dof(i) += epsilon2 * (cos(theta2) * eigenvectors[eigenmode2][i] -
                            sin(theta2) * eigenvectors[conjugate2][i]);
    }
    doc_solution();
    std::cout << "Finished perturbation" << std::endl;
  }

  void solve_for_eigenproblem(Vector<complex<double>>& eigenvalues,
                              Vector<DoubleVector>& eigenvectors)

  {
    bool doc_eigenvector = false;
    bool silent = false;

    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    if (doc_eigenvector == true)
    {
      if (silent == false)
      {
        std::cout << "Doc Initial Solution" << std::endl;
      }

      doc_solution();
    }
    // Reset Eigenvalues

    for (unsigned n = 0; n < 10; n++)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
      Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
    }

    // Vector<complex<double> > eigenvalues;
    // Vector<DoubleVector> eigenvectors;

    unsigned n_eval = 10;
    int sensible_eigenvalues = 0;
    // Solve for eigenproblem
    if (silent == false)
    {
      std::cout << "Now attempt to solve the eigenproblem" << std::endl;
    }

    solve_eigenproblem(n_eval, eigenvalues, eigenvectors);

    if (silent == false)
    {
      std::cout << "N_eval is " << n_eval << std::endl;

      std::cout << "Eigenvalues size is" << eigenvalues.size() << std::endl;

      // Describe Eigenvalues
      std::cout << "And describe eigenvalues" << std::endl;
    }

    n_eval = eigenvalues.size();
    int Number_of_negative_eigenvalues = 0;
    for (unsigned n = 0; n < n_eval; n++)
    {
      if (isinf(real(eigenvalues[n])) != true &&
          isnan(real(eigenvalues[n])) != true)
      {
        if (silent == false)
        {
          std::cout << "Eigenvalue " << eigenvalues[n] << std::endl;
        }
        sensible_eigenvalues++;
        if (real(eigenvalues[n]) < 0)
        {
          Number_of_negative_eigenvalues++;
        }
      }
      if (doc_eigenvector == true)
      {
        if (silent == false)
        {
          std::cout << "Doc eigenvector with eigenvalue" << real(eigenvalues[n])
                    << " " << imag(eigenvalues[n]) << " " << std::endl;
        }

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

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    for (unsigned i = 0; i < n_dof; i++)
    {
      backup_1[i] = dof(i);
    }
    doc_solution();

    for (unsigned m = 0; m < n_eval; m++)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = 100 * eigenvectors[m][i];
      }
      output_eigenvector_boundaries(m, false);
    }


    for (unsigned m = 0; m < n_eval; m++)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = backup_1[i] + 10 * eigenvectors[m][i];
      }
      std::cout << "=====================" << std::endl;
      std::cout << "=====================" << std::endl;
      std::cout << "Perturbing steady soln" << std::endl;
      std::cout << "with eigenvalue " << std::endl;
      std::cout << real(eigenvalues[m]) << " + i" << imag(eigenvalues[m])
                << std::endl;
      std::cout << "=====================" << std::endl;
      std::cout << "=====================" << std::endl;
      output_eigenvector_boundaries(m, true);
    }

    output_eigenvectors_and_eigenvalues(eigenvalues, eigenvectors);
    int eigenmode;
    double epsilon;
    double theta;
    int conjugate;
    bool answer = false;

    for (unsigned m = 0; m < 1000; m++)
    {
      if (answer == false)
      {
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "Choices, choices....." << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "What eigenmode? (0-9)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> eigenmode;
        std::cout << "What is the conjugate? (0-9)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> conjugate;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "What phase? (as a multiple of pi/12)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> theta;
        std::cout << "What perturbation? (any double)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> epsilon;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "Sure you're happy? (true/false)" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << "=====================" << std::endl;
        std::cin >> answer;
      }
      else
      {
        break;
      }
    }

    theta = (1.0 / 12.0) * theta * MathematicalConstants::Pi;

    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "Okay let's do this" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "Perturbing steady soln" << std::endl;
    std::cout << "with eigenvalue " << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << real(eigenvalues[eigenmode]) << " + i"
              << imag(eigenvalues[eigenmode]) << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "with phase" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << theta << std::endl;


    ofstream some_file;
    char filename[100];

    sprintf(filename, "perturbation.dat");

    some_file.open(filename);
    some_file << eigenmode << " " << epsilon << " "
              << " " << theta << endl;
    some_file.close();

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) = 100 * (cos(theta) * eigenvectors[eigenmode][i] -
                      sin(theta) * eigenvectors[conjugate][i]);
    }

    output_eigenvector_boundaries(20, false);

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) =
        backup_1[i] + epsilon * (cos(theta) * eigenvectors[eigenmode][i] -
                                 sin(theta) * eigenvectors[conjugate][i]);
    }
  }

  void solve_for_eigenproblem_edge(double eps)
  {
    Vector<complex<double>> eigenvalues;
    Vector<DoubleVector> eigenvectors;

    solve_for_eigenproblem_edge(eigenvalues, eigenvectors, eps);
  }


  void solve_for_eigenproblem_edge(Vector<complex<double>>& eigenvalues,
                                   Vector<DoubleVector>& eigenvectors,
                                   double eps)

  {
    bool doc_eigenvector = false;
    bool silent = false;

    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    if (doc_eigenvector == true)
    {
      if (silent == false)
      {
        std::cout << "Doc Initial Solution" << std::endl;
      }

      doc_solution();
    }
    // Reset Eigenvalues

    for (unsigned n = 0; n < 10; n++)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
      Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
    }

    // Vector<complex<double> > eigenvalues;
    // Vector<DoubleVector> eigenvectors;

    unsigned n_eval = 10;
    int sensible_eigenvalues = 0;
    // Solve for eigenproblem
    if (silent == false)
    {
      std::cout << "Now attempt to solve the eigenproblem" << std::endl;
    }

    solve_eigenproblem(n_eval, eigenvalues, eigenvectors);

    if (silent == false)
    {
      std::cout << "N_eval is " << n_eval << std::endl;

      std::cout << "Eigenvalues size is" << eigenvalues.size() << std::endl;

      // Describe Eigenvalues
      std::cout << "And describe eigenvalues" << std::endl;
    }

    n_eval = eigenvalues.size();
    int Number_of_negative_eigenvalues = 0;
    for (unsigned n = 0; n < n_eval; n++)
    {
      if (isinf(real(eigenvalues[n])) != true &&
          isnan(real(eigenvalues[n])) != true)
      {
        if (silent == false)
        {
          std::cout << "Eigenvalue " << eigenvalues[n] << std::endl;
        }
        sensible_eigenvalues++;
        if (real(eigenvalues[n]) < 0)
        {
          Number_of_negative_eigenvalues++;
        }
      }
      if (doc_eigenvector == true)
      {
        if (silent == false)
        {
          std::cout << "Doc eigenvector with eigenvalue" << real(eigenvalues[n])
                    << " " << imag(eigenvalues[n]) << " " << std::endl;
        }

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

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    for (unsigned i = 0; i < n_dof; i++)
    {
      backup_1[i] = dof(i);
    }
    doc_solution();

    for (unsigned m = 0; m < n_eval; m++)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = 100 * eigenvectors[m][i];
      }
      output_eigenvector_boundaries(m, false);
    }


    for (unsigned m = 0; m < n_eval; m++)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = backup_1[i] + 10 * eigenvectors[m][i];
      }
      std::cout << "=====================" << std::endl;
      std::cout << "=====================" << std::endl;
      std::cout << "Perturbing steady soln" << std::endl;
      std::cout << "with eigenvalue " << std::endl;
      std::cout << real(eigenvalues[m]) << " + i" << imag(eigenvalues[m])
                << std::endl;
      std::cout << "=====================" << std::endl;
      std::cout << "=====================" << std::endl;
      output_eigenvector_boundaries(m, true);
    }

    output_eigenvectors_and_eigenvalues(eigenvalues, eigenvectors);
    int eigenmode;
    double theta = 1.0;

    theta = (1.0 / 12.0) * theta * MathematicalConstants::Pi;

    ofstream some_file;
    char filename[100];

    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "perturbation.dat").c_str());

    some_file.open(filename);
    some_file << eigenmode << " " << eps << " "
              << " " << theta << endl;
    some_file.close();

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) = backup_1[i] + eps * (cos(theta) * eigenvectors[0][i] -
                                    sin(theta) * eigenvectors[1][i]);
    }
  }


  void solve_for_eigenproblem_fold(Vector<complex<double>>& eigenvalues,
                                   Vector<DoubleVector>& eigenvectors)
  {
    U_data_pt->unpin(0);
    Q_inv_data_pt->pin(0);
    assign_eqn_numbers();


    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "Starting Eigenproblem" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << "=====================" << std::endl;
    std::cout << " " << std::endl;
    bool doc_eigenvector = false;
    bool retain_perturbed_solution_at_function_exit = false;

    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    if (doc_eigenvector == true)
    {
      std::cout << "Doc Initial Solution" << std::endl;
      doc_solution();
    }
    // Reset Eigenvalues

    for (unsigned n = 0; n < 10; n++)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
      Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
    }

    // Vector<complex<double> > eigenvalues;
    // Vector<DoubleVector> eigenvectors;

    unsigned n_eval = 10;
    int sensible_eigenvalues = 0;
    // Solve for eigenproblem

    std::cout << "Now attempt to solve the eigenproblem" << std::endl;
    solve_eigenproblem(n_eval, eigenvalues, eigenvectors);
    std::cout << "N_eval is " << n_eval << std::endl;
    std::cout << "Eigenvalues size is " << eigenvalues.size() << std::endl;

    // Describe Eigenvalues
    std::cout << "And describe eigenvalues" << std::endl;
    n_eval = eigenvalues.size();
    int Number_of_negative_eigenvalues = 0;
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
        std::cout << "Doc eigenvector with eigenvalue" << real(eigenvalues[n])
                  << " " << imag(eigenvalues[n]) << " " << std::endl;

        for (unsigned i = 0; i < n_dof; i++)
        {
          dof(i) = backup_1[i] + 10 * eigenvectors[n][i];
        }
        actions_after_change_in_bifurcation_parameter();
        doc_solution();
      }
    }

    for (unsigned n = 0; n < 10; n++)
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

    std::cout << "There are " << sensible_eigenvalues
              << " reasonable eigenvalues" << std::endl;

    if (retain_perturbed_solution_at_function_exit == true)
    {
      std::cout << "Keep perturbation with eigenvector 0" << std::endl;
      unsigned perturbation_evec = 0;

      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = backup_1[i] + 2 * eigenvectors[perturbation_evec][i];
      }
      actions_after_change_in_bifurcation_parameter();
      // doc_solution();
    };

    U_data_pt->pin(0);
    Q_inv_data_pt->unpin(0);
    assign_eqn_numbers();

    std::cout << "======================" << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Finishing Eigenproblem" << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "======================" << std::endl;
  }

  double solve_for_hopf_bisection()
  {
    // U_data_pt->pin(0);
    // Q_inv_data_pt->unpin(0);
    // set_U(1.9615);
    unsigned eval_index;
    Vector<complex<double>> eigenvalues;
    Vector<DoubleVector> eigenvectors;
    steady_newton_solve();
    solve_for_eigenproblem_fold(eigenvalues, eigenvectors);
    double real_evalue_min = 0;
    double imag_evalue_min = 0;
    eval_index = sort_complex_eigenvalues(&real_evalue_min, &imag_evalue_min);
    double real_evalue_temp = real_evalue_min;
    double imag_evalue_temp = imag_evalue_min;
    std::cout << "===========================" << std::endl;
    std::cout << "Old U is " << get_U() << std::endl;
    std::cout << "Temp Real part is " << real_evalue_temp << std::endl;
    std::cout << "Temp Imag part is " << imag_evalue_temp << std::endl;
    std::cout << "===========================" << std::endl;

    // Vary U

    double mult = 1.0; // If H1 then set mult = 1.0, if H2 then set mult = -1.0;

    double U = get_U();
    // double ds = -0.001;
    double tol = 1e-8;
    double U_old = U;
    double ds = 0.005;
    if (real_evalue_min > 0)
    {
      U = U - mult * ds;
    }
    else
    {
      U = U + mult * ds;
    }

    set_U(U);

    for (unsigned m = 0; m < 1000; m++)
    {
      // Set new value of U

      if (m % 5 == 4)
      {
        steady_newton_solve();
      }
      else
      {
        steady_newton_solve();
      }
      solve_for_eigenproblem_fold(eigenvalues, eigenvectors);
      eval_index = sort_complex_eigenvalues(&real_evalue_min, &imag_evalue_min);
      std::cout << "===========================" << std::endl;
      std::cout << "New U is " << get_U() << std::endl;
      std::cout << "New Real part is " << real_evalue_min << std::endl;
      std::cout << "New Imag part is " << imag_evalue_min << std::endl;
      std::cout << "===========================" << std::endl;
      std::cout << "===========================" << std::endl;
      std::cout << "Old U is " << U_old << std::endl;
      std::cout << "Old Real part is " << real_evalue_temp << std::endl;
      std::cout << "Old Imag part is " << imag_evalue_temp << std::endl;
      std::cout << "===========================" << std::endl;

      // Decide if new value of U is better
      if ((real_evalue_min < 0) & ((real_evalue_min * real_evalue_temp) > 0))
      {
        // If not change step direction of U
        // Return to previous value of U
        U_old = U;
        U = U + mult * ds;
        set_U(U);
        // Update values to better values
        real_evalue_temp = real_evalue_min;
        imag_evalue_temp = imag_evalue_min;
        std::cout << "UNSTABLE" << std::endl;
      }
      if ((real_evalue_min > 0) & ((real_evalue_min * real_evalue_temp) > 0))
      {
        // Otherwise keep the same step
        U_old = U;
        U = U - mult * ds;
        set_U(U);
        // Update values to better values
        real_evalue_temp = real_evalue_min;
        imag_evalue_temp = imag_evalue_min;
        std::cout << "STABLE" << std::endl;
      }

      // Check for a change in sign
      if ((real_evalue_min) * (real_evalue_temp) < 0)
      {
        // If a change in sign then change direction of U step and make smaller
        std::cout << "Change is good for the soul..." << std::endl;
        break;
      }

      // Output Eigenvalue
      reset_lagrangian_coordinates();
    }

    std::cout << "=================" << std::endl;
    std::cout << "Starting Bisection" << std::endl;
    std::cout << "=================" << std::endl;

    // Start Bisection
    double U_neg = 0.0;
    double U_pos = 0.0;
    double U_mid;
    double evalue_pos;
    double evalue_neg;

    // Assign which U is the 'positive'/'negative' one
    if (real_evalue_temp > 0)
    {
      U_pos = U_old;
      U_neg = U;
      evalue_pos = real_evalue_temp;
      evalue_neg = real_evalue_min;
    }

    if (real_evalue_min > 0)
    {
      U_pos = U;
      U_neg = U_old;
      evalue_pos = real_evalue_min;
      evalue_neg = real_evalue_temp;
    }

    // Reset temporary values
    real_evalue_temp = 0.0;
    imag_evalue_temp = 0.0;

    for (unsigned m = 0; m < 100; m++)
    {
      U_mid = 0.5 * (U_neg + U_pos);
      std::cout << "U_mid is " << U_mid << std::endl;

      set_U(U_mid);
      if (m % 5 == 4)
      {
        steady_newton_solve(0);
      }
      else
      {
        steady_newton_solve(0);
      }

      solve_for_eigenproblem_fold(eigenvalues, eigenvectors);
      eval_index =
        sort_complex_eigenvalues(&real_evalue_temp, &imag_evalue_temp);
      reset_lagrangian_coordinates();
      // Decide if the imaginary part is zero
      if (sqrt(real_evalue_temp * real_evalue_temp) < tol)
      {
        break;
      }

      // If not update the values of U

      if (real_evalue_temp > 0)
      {
        evalue_pos = real_evalue_temp;
        U_pos = U_mid;
      }

      if (real_evalue_temp < 0)
      {
        evalue_neg = real_evalue_temp;
        U_neg = U_mid;
      }

      // If U is correct to 7 s.f
      if (sqrt((U_pos - U_neg) * (U_pos - U_neg)) < tol)
      {
        break;
      }

      std::cout << "=======================" << std::endl;
      std::cout << "U_pos is " << U_pos << std::endl;
      std::cout << "evalue_pos is " << evalue_pos << std::endl;
      std::cout << "U_neg is " << U_neg << std::endl;
      std::cout << "evalue_neg is " << evalue_neg << std::endl;
      std::cout << "=======================" << std::endl;
    }


    // Bifurcation at
    // doc_solution();
    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    // Doc the eigenvectors

    std::cout << "Doc eigenvector with eigenvalue"
              << real(eigenvalues[eval_index]) << " "
              << imag(eigenvalues[eval_index]) << " " << std::endl;

    /* for(unsigned i=0;i<n_dof;i++) */
    /* 	{ */
    /* 	  dof(i) = backup_1[i]+10* eigenvectors[eval_index][i]; */

    /* 	} */
    /* doc_solution(); */
    /* for(unsigned i=0;i<n_dof;i++) */
    /* 	{ */
    /* 	  dof(i) = backup_1[i]; */
    /* 	} */
    doc_solution();
    output_eigenvectors_and_eigenvalues(eigenvalues, eigenvectors);

    for (unsigned i = 0; i < n_dof; i++)
    {
      backup_1[i] = dof(i);
    }

    double n_eval = eigenvalues.size();

    for (unsigned m = 0; m < n_eval; m++)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = 100 * eigenvectors[m][i];
      }
      output_eigenvector_boundaries(m);
    }

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) = backup_1[i];
    }

    std::cout << "========================" << std::endl;
    std::cout << "==BIFURCATION LOCATION==" << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << "Q is " << get_Q() << std::endl;
    std::cout << "U is " << get_U() << std::endl;
    std::cout << "V is " << get_V() << std::endl;
    std::cout << "Eigenvalue is " << real_evalue_temp << " + i"
              << imag_evalue_temp << std::endl;
    std::cout << "================" << std::endl;
    Q_inv_data_pt->pin(0);
    U_data_pt->unpin(0);
    assign_eqn_numbers();

    return ds;
  }

  double solve_for_fold_bisection()
  {
    // U_data_pt->pin(0);
    // Q_inv_data_pt->unpin(0);
    // set_U(1.9615);
    unsigned eval_index;
    Vector<complex<double>> eigenvalues;
    Vector<DoubleVector> eigenvectors;
    steady_newton_solve();
    solve_for_eigenproblem_fold(eigenvalues, eigenvectors);
    double real_evalue_min = 0;
    double imag_evalue_min = 0;
    eval_index = sort_real_eigenvalues(&real_evalue_min, &imag_evalue_min);
    double real_evalue_temp = real_evalue_min;
    double imag_evalue_temp = imag_evalue_min;
    std::cout << "===========================" << std::endl;
    std::cout << "Old U is " << get_U() << std::endl;
    std::cout << "Temp Real part is " << real_evalue_temp << std::endl;
    std::cout << "Temp Imag part is " << imag_evalue_temp << std::endl;
    std::cout << "===========================" << std::endl;

    // Vary U

    double mult = 1.0; // If H1 then set mult = 1.0, if H2 then set mult = -1.0;

    double U = get_U();
    // double ds = -0.001;
    double tol = 1e-8;
    double U_old = U;
    double ds = 0.005;
    if (real_evalue_min > 0)
    {
      U = U - mult * ds;
    }
    else
    {
      U = U + mult * ds;
    }

    set_U(U);

    for (unsigned m = 0; m < 1000; m++)
    {
      // Set new value of U

      if (m % 5 == 4)
      {
        steady_newton_solve();
      }
      else
      {
        steady_newton_solve();
      }
      solve_for_eigenproblem_fold(eigenvalues, eigenvectors);
      eval_index = sort_real_eigenvalues(&real_evalue_min, &imag_evalue_min);
      std::cout << "===========================" << std::endl;
      std::cout << "New U is " << get_U() << std::endl;
      std::cout << "New Real part is " << real_evalue_min << std::endl;
      std::cout << "New Imag part is " << imag_evalue_min << std::endl;
      std::cout << "===========================" << std::endl;
      std::cout << "===========================" << std::endl;
      std::cout << "Old U is " << U_old << std::endl;
      std::cout << "Old Real part is " << real_evalue_temp << std::endl;
      std::cout << "Old Imag part is " << imag_evalue_temp << std::endl;
      std::cout << "===========================" << std::endl;

      // Decide if new value of U is better
      if ((real_evalue_min < 0) & ((real_evalue_min * real_evalue_temp) > 0))
      {
        // If not change step direction of U
        // Return to previous value of U
        U_old = U;
        U = U + mult * ds;
        set_U(U);
        // Update values to better values
        real_evalue_temp = real_evalue_min;
        imag_evalue_temp = imag_evalue_min;
        std::cout << "UNSTABLE" << std::endl;
      }
      if ((real_evalue_min > 0) & ((real_evalue_min * real_evalue_temp) > 0))
      {
        // Otherwise keep the same step
        U_old = U;
        U = U - mult * ds;
        set_U(U);
        // Update values to better values
        real_evalue_temp = real_evalue_min;
        imag_evalue_temp = imag_evalue_min;
        std::cout << "STABLE" << std::endl;
      }

      // Check for a change in sign
      if ((real_evalue_min) * (real_evalue_temp) < 0)
      {
        // If a change in sign then change direction of U step and make smaller
        std::cout << "Change is good for the soul..." << std::endl;
        break;
      }

      // Output Eigenvalue
      reset_lagrangian_coordinates();
    }

    std::cout << "=================" << std::endl;
    std::cout << "Starting Bisection" << std::endl;
    std::cout << "=================" << std::endl;

    // Start Bisection
    double U_neg = 0.0;
    double U_pos = 0.0;
    double U_mid;
    double evalue_pos;
    double evalue_neg;

    // Assign which U is the 'positive'/'negative' one
    if (real_evalue_temp > 0)
    {
      U_pos = U_old;
      U_neg = U;
      evalue_pos = real_evalue_temp;
      evalue_neg = real_evalue_min;
    }

    if (real_evalue_min > 0)
    {
      U_pos = U;
      U_neg = U_old;
      evalue_pos = real_evalue_min;
      evalue_neg = real_evalue_temp;
    }

    // Reset temporary values
    real_evalue_temp = 0.0;
    imag_evalue_temp = 0.0;

    for (unsigned m = 0; m < 100; m++)
    {
      U_mid = 0.5 * (U_neg + U_pos);
      std::cout << "U_mid is " << U_mid << std::endl;

      set_U(U_mid);
      if (m % 5 == 4)
      {
        steady_newton_solve(0);
      }
      else
      {
        steady_newton_solve(0);
      }

      solve_for_eigenproblem_fold(eigenvalues, eigenvectors);
      eval_index = sort_real_eigenvalues(&real_evalue_temp, &imag_evalue_temp);
      reset_lagrangian_coordinates();
      // Decide if the imaginary part is zero
      if (sqrt(real_evalue_temp * real_evalue_temp) < tol)
      {
        break;
      }

      // If not update the values of U

      if (real_evalue_temp > 0)
      {
        evalue_pos = real_evalue_temp;
        U_pos = U_mid;
      }

      if (real_evalue_temp < 0)
      {
        evalue_neg = real_evalue_temp;
        U_neg = U_mid;
      }

      // If U is correct to 7 s.f
      if (sqrt((U_pos - U_neg) * (U_pos - U_neg)) < tol)
      {
        break;
      }

      std::cout << "=======================" << std::endl;
      std::cout << "U_pos is " << U_pos << std::endl;
      std::cout << "evalue_pos is " << evalue_pos << std::endl;
      std::cout << "U_neg is " << U_neg << std::endl;
      std::cout << "evalue_neg is " << evalue_neg << std::endl;
      std::cout << "=======================" << std::endl;
    }


    // Bifurcation at
    // doc_solution();
    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    // Doc the eigenvectors

    std::cout << "Doc eigenvector with eigenvalue"
              << real(eigenvalues[eval_index]) << " "
              << imag(eigenvalues[eval_index]) << " " << std::endl;

    /* for(unsigned i=0;i<n_dof;i++) */
    /* 	{ */
    /* 	  dof(i) = backup_1[i]+10* eigenvectors[eval_index][i]; */

    /* 	} */
    /* doc_solution(); */
    /* for(unsigned i=0;i<n_dof;i++) */
    /* 	{ */
    /* 	  dof(i) = backup_1[i]; */
    /* 	} */
    doc_solution();
    output_eigenvectors_and_eigenvalues(eigenvalues, eigenvectors);

    for (unsigned i = 0; i < n_dof; i++)
    {
      backup_1[i] = dof(i);
    }

    double n_eval = eigenvalues.size();

    for (unsigned m = 0; m < n_eval; m++)
    {
      for (unsigned i = 0; i < n_dof; i++)
      {
        dof(i) = 100 * eigenvectors[m][i];
      }
      output_eigenvector_boundaries(m);
    }

    for (unsigned i = 0; i < n_dof; i++)
    {
      dof(i) = backup_1[i];
    }

    std::cout << "========================" << std::endl;
    std::cout << "==BIFURCATION LOCATION==" << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << "Q is " << get_Q() << std::endl;
    std::cout << "U is " << get_U() << std::endl;
    std::cout << "V is " << get_V() << std::endl;
    std::cout << "Eigenvalue is " << real_evalue_temp << " + i"
              << imag_evalue_temp << std::endl;
    std::cout << "================" << std::endl;
    return ds;
    steady_newton_solve(1);
  }


  void solve_for_eigenproblem_no_perturb()
  {
    Vector<complex<double>> eigenvalues;
    Vector<DoubleVector> eigenvectors;

    solve_for_eigenproblem_no_perturb(eigenvalues, eigenvectors);
  }


  void solve_for_eigenproblem_no_perturb(Vector<complex<double>>& eigenvalues,
                                         Vector<DoubleVector>& eigenvectors)
  {
    bool doc_eigenvector = false;
    bool silent = false;

    //      Q_inv_data_pt->pin(0);
    // U_data_pt->unpin(0);
    // this->assign_eqn_numbers();

    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    if (doc_eigenvector == true)
    {
      if (silent == false)
      {
        std::cout << "Doc Initial Solution" << std::endl;
      }

      doc_solution();
    }
    // Reset Eigenvalues

    for (unsigned n = 0; n < 10; n++)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
      Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
    }

    // Vector<complex<double> > eigenvalues;
    // Vector<DoubleVector> eigenvectors;

    unsigned n_eval = 10;
    int sensible_eigenvalues = 0;
    // Solve for eigenproblem
    if (silent == false)

    {
      std::cout << "Now attempt to solve the eigenproblem" << std::endl;
    }
    solve_eigenproblem(n_eval, eigenvalues, eigenvectors);

    if (silent == false)
    {
      std::cout << "N_eval is " << n_eval << std::endl;

      std::cout << "Eigenvalues size is" << eigenvalues.size() << std::endl;

      // Describe Eigenvalues
      std::cout << "And describe eigenvalues" << std::endl;
    }

    n_eval = eigenvalues.size();
    int Number_of_negative_eigenvalues = 0;
    for (unsigned n = 0; n < n_eval; n++)
    {
      if (isinf(real(eigenvalues[n])) != true &&
          isnan(real(eigenvalues[n])) != true)
      {
        if (silent == false)
        {
          std::cout << "Eigenvalue " << eigenvalues[n] << std::endl;
        }
        sensible_eigenvalues++;
        if (real(eigenvalues[n]) < 0)
        {
          Number_of_negative_eigenvalues++;
        }
      }
      if (doc_eigenvector == true)
      {
        if (silent == false)
        {
          std::cout << "Doc eigenvector with eigenvalue" << real(eigenvalues[n])
                    << " " << imag(eigenvalues[n]) << " " << std::endl;
        }

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

    output_eigenvectors_and_eigenvalues(eigenvalues, eigenvectors);

    //      Q_inv_data_pt->unpin(0);
    // U_data_pt->pin(0);
    // this->assign_eqn_numbers();
  }

  void solve_for_adjoint_eigenproblem(Vector<complex<double>>& eigenvalues,
                                      Vector<DoubleVector>& eigenvectors)
  {
    bool doc_eigenvector = false;
    bool silent = false;

    //      Q_inv_data_pt->pin(0);
    // U_data_pt->unpin(0);
    // this->assign_eqn_numbers();

    unsigned n_dof = ndof();
    Vector<double> backup_1(n_dof);

    // Back up the solution
    for (unsigned n = 0; n < n_dof; n++)
    {
      backup_1[n] = dof(n);
    }

    if (doc_eigenvector == true)
    {
      if (silent == false)
      {
        std::cout << "Doc Initial Solution" << std::endl;
      }

      doc_solution();
    }
    // Reset Eigenvalues

    for (unsigned n = 0; n < 10; n++)
    {
      Problem_Parameter::vector_of_eigenvalues_rp[n] = 10;
      Problem_Parameter::vector_of_eigenvalues_ip[n] = 10;
    }

    // Vector<complex<double> > eigenvalues;
    // Vector<DoubleVector> eigenvectors;

    unsigned n_eval = 10;
    int sensible_eigenvalues = 0;
    // Solve for eigenproblem
    if (silent == false)

    {
      std::cout << "Now attempt to solve the eigenproblem" << std::endl;
    }
    solve_adjoint_eigenproblem(n_eval, eigenvalues, eigenvectors);

    if (silent == false)
    {
      std::cout << "N_eval is " << n_eval << std::endl;

      std::cout << "Eigenvalues size is" << eigenvalues.size() << std::endl;

      // Describe Eigenvalues
      std::cout << "And describe eigenvalues" << std::endl;
    }

    n_eval = eigenvalues.size();
    int Number_of_negative_eigenvalues = 0;
    for (unsigned n = 0; n < n_eval; n++)
    {
      if (isinf(real(eigenvalues[n])) != true &&
          isnan(real(eigenvalues[n])) != true)
      {
        if (silent == false)
        {
          std::cout << "Eigenvalue " << eigenvalues[n] << std::endl;
        }
        sensible_eigenvalues++;
        if (real(eigenvalues[n]) < 0)
        {
          Number_of_negative_eigenvalues++;
        }
      }
      if (doc_eigenvector == true)
      {
        if (silent == false)
        {
          std::cout << "Doc eigenvector with eigenvalue" << real(eigenvalues[n])
                    << " " << imag(eigenvalues[n]) << " " << std::endl;
        }

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

    //      Q_inv_data_pt->unpin(0);
    // U_data_pt->pin(0);
    // this->assign_eqn_numbers();
  }

  void output_eigenvector_perturbation(int eps, int theta)
  {
    ofstream some_file;
    char filename[100];
    sprintf(filename,
            "%s/perturbation_eps_%i_theta_%i.dat",
            Problem_Parameter::Doc_info.directory().c_str(),
            eps,
            theta);

    some_file.open(filename);
    this->Fluid_mesh_pt->output_boundaries(some_file);
    some_file.close();
  }

  void output_eigenvector_boundaries(unsigned m, bool pert = false)
  {
    ofstream some_file;
    char filename[100];

    if (pert == false)
    {
      sprintf(filename,
              "%s/boundaries_evectors_%i_%i.dat",
              Problem_Parameter::Doc_info.directory().c_str(),
              Problem_Parameter::Doc_info.number(),
              m);
    }

    if (pert == true)
    {
      sprintf(filename,
              "%s/boundaries_pert_evectors_%i_%i.dat",
              Problem_Parameter::Doc_info.directory().c_str(),
              Problem_Parameter::Doc_info.number(),
              m);
    }

    some_file.open(filename);
    this->Fluid_mesh_pt->output_boundaries(some_file);
    some_file.close();
  }

  void output_eigenvectors_and_eigenvalues(Vector<complex<double>>& eigenvalues,
                                           Vector<DoubleVector>& eigenvectors)
  {
    double n_evals = eigenvalues.size();
    double n_dof = ndof();
    ofstream some_file;
    char filename[100];

    for (int i = 0; i < n_evals; i++)
    {
      sprintf(filename,
              (Problem_Parameter::Doc_info.directory() + "eigenvectors_%i.dat")
                .c_str(),
              i);
      some_file.open(filename);

      for (int j = 0; j < n_dof; j++)
      {
        some_file << eigenvectors[i][j] << " ";
      }

      some_file << std::endl;
      some_file.close();
    }


    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvalues_%i.dat").c_str(),
      Problem_Parameter::Doc_info.number());
    some_file.open(filename);

    for (int i = 0; i < n_evals; i++)
    {
      some_file << real(eigenvalues[i]) << " " << imag(eigenvalues[i]) << " ";
      some_file << std::endl;
    }

    some_file.close();
  }

  void output_adjoint_eigenvectors_and_eigenvalues(
    Vector<complex<double>>& eigenvalues, Vector<DoubleVector>& eigenvectors)
  {
    double n_evals = eigenvalues.size();
    double n_dof = ndof();
    ofstream some_file;
    char filename[100];

    for (int i = 0; i < n_evals; i++)
    {
      sprintf(filename,
              (Problem_Parameter::Doc_info.directory() +
               "eigenvectors_adjoint_%i.dat")
                .c_str(),
              i);
      some_file.open(filename);

      for (int j = 0; j < n_dof; j++)
      {
        some_file << eigenvectors[i][j] << " ";
      }

      some_file << std::endl;
      some_file.close();
    }


    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "eigenvalues_adjoint_%i.dat")
        .c_str(),
      Problem_Parameter::Doc_info.number());
    some_file.open(filename);

    for (int i = 0; i < n_evals; i++)
    {
      some_file << real(eigenvalues[i]) << " " << imag(eigenvalues[i]) << " ";
      some_file << std::endl;
    }

    some_file.close();
  }


  void output_eigenmatrices()
  {
    std::cout << "OUTPUT EIGENMATRICES" << std::endl;
    CRDoubleMatrix temp_J(this->dof_distribution_pt());
    CRDoubleMatrix temp_M(this->dof_distribution_pt());

    std::cout << "Get mass matrix " << std::endl;
    Problem::get_eigenproblem_matrices(temp_M, temp_J);

    int n_dof = ndof();
    ofstream some_file;
    char filename[100];

    std::cout << "Now write mass matrix to file " << std::endl;
    sprintf(
      filename,
      (Problem_Parameter::Doc_info.directory() + "mass_matrix_%i.dat").c_str(),
      Problem_Parameter::Doc_info.number());

    some_file.open(filename);
    for (int i = 0; i < n_dof; i++)
    {
      for (int j = 0; j < n_dof; j++)
      {
        some_file << temp_M(i, j) << " ";
      }
      some_file << std::endl;
    }

    std::cout << "Now write jacobian matrix to file " << std::endl;
    sprintf(filename,
            (Problem_Parameter::Doc_info.directory() + "jacobian_matrix_%i.dat")
              .c_str(),
            Problem_Parameter::Doc_info.number());

    some_file.open(filename);
    for (int i = 0; i < n_dof; i++)
    {
      for (int j = 0; j < n_dof; j++)
      {
        some_file << temp_J(i, j) << " ";
      }
      some_file << std::endl;
    }
  }

  unsigned sort_eigenvalues(double* real_evalue_min, double* imag_evalue_min)
  {
    // Finds the eigenvalue with minimum real part.

    double real_evalue_temp = Problem_Parameter::vector_of_eigenvalues_rp[0];
    unsigned evalue_index = 0;

    for (unsigned m = 1; m < 10; m++)
    {
      if (Problem_Parameter::vector_of_eigenvalues_rp[m] < real_evalue_temp)
      {
        real_evalue_temp = Problem_Parameter::vector_of_eigenvalues_rp[m];
        evalue_index = m;
      }
      else
      {
      }
    }

    *imag_evalue_min =
      Problem_Parameter::vector_of_eigenvalues_ip[evalue_index];
    *real_evalue_min =
      Problem_Parameter::vector_of_eigenvalues_rp[evalue_index];
    std::cout << "===========================" << std::endl;
    std::cout << "Smallest eigenvalue is at " << evalue_index << std::endl;
    std::cout << "Real part is " << *real_evalue_min << std::endl;
    std::cout << "Imag part is " << *imag_evalue_min << std::endl;
    std::cout << "===========================" << std::endl;
    return evalue_index;
  }


  unsigned sort_real_eigenvalues(double* real_evalue_min,
                                 double* imag_evalue_min)
  {
    // Finds the eigenvalue with minimum real part but only real eigenvalues

    int real_index = 0;
    int counter = 0;

    for (int m = 0; m < 10; m++)
    {
      if (pow(Problem_Parameter::vector_of_eigenvalues_ip[m], 2) < 1e-9)
      {
        real_index += 1;
      }
      else
      {
      }
    }

    Vector<int> real(real_index);

    for (int m = 0; m < 10; m++)
    {
      if (pow(Problem_Parameter::vector_of_eigenvalues_ip[m], 2) < 1e-9)
      {
        real[counter] = m;
        counter += 1;
      }
    }

    std::cout << " " << std::endl;
    std::cout << "There are " << real_index << " eigenvalues which are real"
              << std::endl;
    std::cout << "They are at.." << std::endl;
    for (int m = 0; m < real_index; m++)
    {
      std::cout << real[m] << std::endl;
    }
    std::cout << " " << std::endl;

    double real_evalue_temp =
      Problem_Parameter::vector_of_eigenvalues_rp[real[0]];
    unsigned evalue_index = real[0];

    for (int m = 1; m < real_index; m++)
    {
      if (Problem_Parameter::vector_of_eigenvalues_rp[real[m]] <
          real_evalue_temp)
      {
        real_evalue_temp = Problem_Parameter::vector_of_eigenvalues_rp[real[m]];
        evalue_index = real[m];
      }
      else
      {
      }
    }

    *imag_evalue_min =
      Problem_Parameter::vector_of_eigenvalues_ip[evalue_index];
    *real_evalue_min =
      Problem_Parameter::vector_of_eigenvalues_rp[evalue_index];
    std::cout << "===========================" << std::endl;
    std::cout << "Smallest eigenvalue is at " << evalue_index << std::endl;
    std::cout << "Real part is " << *real_evalue_min << std::endl;
    std::cout << "Imag part is " << *imag_evalue_min << std::endl;
    std::cout << "===========================" << std::endl;
    return evalue_index;
  }

  unsigned sort_complex_eigenvalues(double* real_evalue_min,
                                    double* imag_evalue_min)
  {
    // Finds the eigenvalue with minimum real part but only real eigenvalues

    int complex_index = 0;
    int counter = 0;

    for (int m = 0; m < 10; m++)
    {
      if (Problem_Parameter::vector_of_eigenvalues_ip[m] != 0)
      {
        complex_index += 1;
      }
      else
      {
      }
    }

    Vector<int> imag(complex_index);

    for (int m = 0; m < 10; m++)
    {
      if (Problem_Parameter::vector_of_eigenvalues_ip[m] != 0)
      {
        imag[counter] = m;
        counter += 1;
      }
    }

    std::cout << " " << std::endl;
    std::cout << "There are " << complex_index
              << " eigenvalues with an imaginary part" << std::endl;
    std::cout << "They are at.." << std::endl;
    for (int m = 0; m < complex_index; m++)
    {
      std::cout << imag[m] << std::endl;
    }
    std::cout << " " << std::endl;


    double complex_evalue_temp =
      Problem_Parameter::vector_of_eigenvalues_rp[imag[0]];
    unsigned evalue_index = imag[0];

    for (unsigned m = 1; m < complex_index; m++)
    {
      if (Problem_Parameter::vector_of_eigenvalues_rp[imag[m]] <
          complex_evalue_temp)
      {
        complex_evalue_temp =
          Problem_Parameter::vector_of_eigenvalues_rp[imag[m]];
        evalue_index = imag[m];
      }
      else
      {
      }
    }

    *imag_evalue_min =
      Problem_Parameter::vector_of_eigenvalues_ip[evalue_index];
    *real_evalue_min =
      Problem_Parameter::vector_of_eigenvalues_rp[evalue_index];
    std::cout << "===========================" << std::endl;
    std::cout << "Smallest eigenvalue is at " << evalue_index << std::endl;
    std::cout << "Real part is " << *real_evalue_min << std::endl;
    std::cout << "Imag part is " << *imag_evalue_min << std::endl;
    std::cout << "===========================" << std::endl;
    return evalue_index;
  }


  /// \short Set boundary conditions and complete the build of all elements
  void complete_problem_setup();

  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  void doc_solution_trace(const std::string& comment = "");

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
      el_pt->flux_fct_pt() =
        &Problem_Parameter::
          normal_flux_behind_bubble; // JACK - THIS SAYS THAT THE flux_fc_pt() =
                                     // normal_flux_behind_bubble

      /// This one is important!
      el_pt->add_external_data(
        G_data_pt,
        true); // JACK - THIS TELLS THE FLUX ELEMENT ABOUT THE PRESCRIBED G
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

  /// Imposed flow rate
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

  double get_G()
  {
    return G_data_pt->value(0);
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

  void create_CoM_X_constraint_for_Q_elements();
  void delete_CoM_X_constraint_for_Q_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = CoM_X_constraint_for_Q_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 1; e < n_element; e++)
    {
      // Kill surface element
      delete CoM_X_constraint_for_Q_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    CoM_X_constraint_for_Q_mesh_pt->flush_element_and_node_storage();
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
  Mesh* CoM_X_constraint_for_Q_mesh_pt;
  Mesh* CoM_Y_constraint_mesh_pt;
  Mesh* Integral_measures_mesh_pt;

  /// Pointer to Fluid_mesh
  RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;


  // Jack - this is where the Bubble_polygon_pt is defined and it is public to
  // the problem

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
  Data* Q_zero_data_pt;

  ///---------------------
  Data* Asymmetry_data_pt;
  Data* CoM_Y_data_pt;

  Data* Integral_measures_data_pt;

  Mesh* Inflow_mesh_pt;
  Mesh* Outflow_mesh_pt;
  //// /// Pointer to element that imposes volume constraint for bubble
  VolumeConstraintElement* Vol_constraint_el_pt;
  VolumeConstraintElement* CoM_X_constraint_el_pt;
  VolumeConstraintElement* CoM_X_constraint_for_Q_el_pt;
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
  Desired_newton_iterations_ds = 1;
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
  CoM_Y_data_pt->set_value(0, 0.0);
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

  // JACK - The next loop makes the integral measures free parameters. Pin one
  // of them if a particular branch is required, i.e. let CoM_Y = 0 to get a
  // branch that is symmetric - note that this can also be achieved above

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
  // JACK - NEED TO FIND OUT WHAT VOLUME CONSTRAINT ELEMENT DOES - IT MUST ADD
  // VOLUME TO THE RESIDUAL
  CoM_X_constraint_el_pt = new VolumeConstraintElement(
    &Problem_Parameter::Centre_of_mass, U_data_pt, index_of_traded_pressure);
  CoM_X_constraint_for_Q_el_pt =
    new VolumeConstraintElement(&Problem_Parameter::Centre_of_mass,
                                Q_inv_data_pt,
                                index_of_traded_pressure);

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

  // Jack now need to say that the internal point isn't fixed


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
  Fluid_mesh_pt->max_permitted_error() = 2e-5; /// 5e-2;2e-5
  Fluid_mesh_pt->min_permitted_error() = 5e-6; /// 1e-4;5e-6
  Fluid_mesh_pt->max_element_size() = 0.05; /// 0.2;0.05
  Fluid_mesh_pt->min_element_size() = 1e-6; /// 4e-5;1e-6

  // Set boundary condition and complete the build of all elements
  complete_problem_setup();

  // Construct the mesh of free surface elements
  Free_surface_mesh_pt = new Mesh;
  create_free_surface_elements();

  Volume_constraint_mesh_pt = new Mesh;
  create_volume_constraint_elements();


  CoM_X_constraint_mesh_pt = new Mesh;
  create_CoM_X_constraint_elements();

  CoM_X_constraint_for_Q_mesh_pt = new Mesh;
  create_CoM_X_constraint_for_Q_elements();

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
  this->add_sub_mesh(this->CoM_X_constraint_for_Q_mesh_pt);
  this->add_sub_mesh(this->CoM_Y_constraint_mesh_pt);
  this->add_sub_mesh(this->Integral_measures_mesh_pt);

  this->add_sub_mesh(this->Inflow_mesh_pt);
  //    this->add_sub_mesh(this->Outflow_mesh_pt);

  // Build global mesh
  this->build_global_mesh();

  pin_tangential_lagrange();

  Use_finite_differences_for_continuation_derivatives = true;

  Max_residuals = 3000000;
  Max_newton_iterations = 100;
  // Setup equation numbering scheme
  cout << "Number of equations: " << this->assign_eqn_numbers() << std::endl;
//  linear_solver_pt()=new FD_LU;
// Choose eigensolver
#ifdef OOMPH_HAS_TRILINOS
  eigen_solver_pt() = new ANASAZI;
  static_cast<ANASAZI*>(eigen_solver_pt())->set_shift(-10);
#endif

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

      // Jack - it is in the following lines that the hele_shaw_interface_file
      // interacts with all of the parameters and functions provided in the
      // Problem_Parameter namespace file etc.


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


template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_CoM_X_constraint_for_Q_elements()
{
  // Add volume constraint element to the mesh
  CoM_X_constraint_for_Q_mesh_pt->add_element_pt(CoM_X_constraint_for_Q_el_pt);


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
      el_pt->set_volume_constraint_element(CoM_X_constraint_for_Q_el_pt);

      // Add it to the mesh
      CoM_X_constraint_for_Q_mesh_pt->add_element_pt(el_pt);
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
void BubbleInChannelProblem<ELEMENT>::jack_solve(double Q_unsteady,
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
  double new_Q = Q_unsteady;
  double new_Q_inv = 1 / Q_unsteady;
  set_Q(new_Q);
  double ds = 0.1;
  for (unsigned m = 0; m < n_iter; m++)
  {
    Q_inv_data_pt->pin(0);

    //      assign_eqn_numbers();
    std::cout << "===========================================" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "CONTINUATION STEP " << m << std::endl;
    std::cout << "Ca is " << get_Q() * get_U() << std::endl;
    std::cout << "U is " << get_U() << std::endl;
    std::cout << "COM_Y is " << get_CoM_Y() << std::endl;
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
        //		  ds = 0.1;
      }
      else
      {
        ds = 0.001;
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

    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "COM_Y is " << get_CoM_Y() << std::endl;
    std::cout << "U_b is " << get_U() << std::endl;
    std::cout << "Ca is " << get_U() * get_Q() << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    doc_solution();
    reset_lagrangian_coordinates();
  }
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::fill_in_landau_coefficients(
  Vector<complex<double>>& g,
  Vector<complex<double>>& g_adjoint,
  Vector<complex<double>>& varphi_0,
  Vector<complex<double>>& varphi_2,
  Vector<complex<double>>& varphi_3,
  double Q_c,
  double omega_c)
{
  int N = this->ndof();
  CRDoubleMatrix jacobian(this->dof_distribution_pt());
  CRDoubleMatrix mass(this->dof_distribution_pt());
  DoubleVector residuals(this->dof_distribution_pt());

  // Lets get the eigenproblem matrices
  Problem::get_eigenproblem_matrices(mass, jacobian);

  // Store the solution
  Vector<double> u0(N);
  for (int i = 0; i < N; i++)
  {
    u0[i] = this->dof(i);
  }

  // Create all of the necessary vectors
  Vector<complex<double>> f(N);
  DoubleVector f_real(this->dof_distribution_pt());
  DoubleVector f_imag(this->dof_distribution_pt());
  DoubleVector g_real(this->dof_distribution_pt());
  DoubleVector g_imag(this->dof_distribution_pt());
  DoubleVector g_conj_real(this->dof_distribution_pt());
  DoubleVector g_conj_imag(this->dof_distribution_pt());
  DoubleVector g_adjoint_real(this->dof_distribution_pt());
  DoubleVector g_adjoint_imag(this->dof_distribution_pt());
  DoubleVector varphi_0_real(this->dof_distribution_pt());
  DoubleVector varphi_0_imag(this->dof_distribution_pt());
  DoubleVector varphi_2_real(this->dof_distribution_pt());
  DoubleVector varphi_2_imag(this->dof_distribution_pt());
  DoubleVector varphi_3_real(this->dof_distribution_pt());
  DoubleVector varphi_3_imag(this->dof_distribution_pt());


  // Initialise the vectors
  for (int i = 0; i < N; i++)
  {
    g_real[i] = real(g[i]);
    g_imag[i] = imag(g[i]);
    g_conj_real[i] = real(g[i]);
    g_conj_imag[i] = -1.0 * imag(g[i]);
    g_adjoint_real[i] = real(g_adjoint[i]);
    g_adjoint_imag[i] = imag(g_adjoint[i]);
    f_real[i] = 0.0;
    f_imag[i] = 0.0;
    varphi_0_real[i] = real(varphi_0[i]);
    varphi_0_imag[i] = imag(varphi_0[i]);
    varphi_2_real[i] = real(varphi_2[i]);
    varphi_2_imag[i] = imag(varphi_2[i]);
    varphi_3_real[i] = real(varphi_3[i]);
    varphi_3_imag[i] = imag(varphi_3[i]);
  }

  // Landau coefficients
  complex<double> lambda = 0.0;
  complex<double> mu = 0.0;
  complex<double> nu = 0.0;

  //================================================================
  // Get the nu coefficient

  // Create complex mass
  CRComplexMatrix mass_complex;
  create_complex_mass_matrix(mass_complex, 1.0);

  // M*g
  mass_complex.multiply(g, f);

  //<M*g,g_adj>
  for (int i = 0; i < N; i++)
  {
    nu += f[i] * g_adjoint[i];
  }

  std::cout << " " << std::endl;
  std::cout << "nu has been calculated " << std::endl;
  //================================================================
  // Get the lambda coefficient

  // Get the Q_c*J[u0]g term

  set_Q(1e30);
  CRDoubleMatrix jacobian_plus(this->dof_distribution_pt());
  Problem::get_jacobian(residuals, jacobian_plus);
  mult_constant(jacobian_plus, -1.0);
  jacobian.add(jacobian_plus, jacobian);

  // Real part
  jacobian.multiply(g_real, f_real);

  // Imag part
  jacobian.multiply(g_imag, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = f_real[i] + 1i * f_imag[i];
  }

  // Update lambda
  for (int i = 0; i < N; i++)
  {
    lambda += f[i] * g_adjoint[i];
  }

  // Reset Q
  set_Q(Q_c);
  Problem::get_jacobian(residuals, jacobian);
  Problem::get_eigenproblem_matrices(mass, jacobian);

  // Get the H[u0](g,varphi_3) term

  get_hessian(g_real, varphi_3_real, u0, f_real);
  get_hessian(g_real, varphi_3_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (f_real[i] + 1i * f_imag[i]);
  }

  // Update lambda
  for (int i = 0; i < N; i++)
  {
    lambda += f[i] * g_adjoint[i];
  }

  get_hessian(g_imag, varphi_3_imag, u0, f_real);
  get_hessian(g_imag, varphi_3_real, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (-1.0 * f_real[i] + 1i * f_imag[i]);
  }

  // Update lambda
  for (int i = 0; i < N; i++)
  {
    lambda += f[i] * g_adjoint[i];
  }

  // Get the H[u0](varphi_3,g) term

  get_hessian(varphi_3_real, g_real, u0, f_real);
  get_hessian(varphi_3_real, g_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (f_real[i] + 1i * f_imag[i]);
  }

  // Update lambda
  for (int i = 0; i < N; i++)
  {
    lambda += f[i] * g_adjoint[i];
  }

  get_hessian(varphi_3_imag, g_real, u0, f_imag);
  get_hessian(varphi_3_imag, g_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (-1.0 * f_real[i] + 1i * f_imag[i]);
  }

  // Update lambda
  for (int i = 0; i < N; i++)
  {
    lambda += f[i] * g_adjoint[i];
  }

  // Get i*omega_c*M[varphu_3]g term

  get_mass(varphi_3_real, g_real, u0, f_real);
  get_mass(varphi_3_real, g_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 1i * omega_c * (f_real[i] + f_imag[i]);
  }

  // Update lambda
  for (int i = 0; i < N; i++)
  {
    lambda += f[i] * g_adjoint[i];
  }

  get_mass(varphi_3_imag, g_real, u0, f_imag);
  get_mass(varphi_3_imag, g_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = 1i * omega_c * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update lambda
  for (int i = 0; i < N; i++)
  {
    lambda += f[i] * g_adjoint[i];
  }

  // Reset f
  for (int i = 0; i < N; i++)
  {
    f[i] = 0.0;
    f_real[i] = 0.0;
    f_imag[i] = 0.0;
  }

  std::cout << " " << std::endl;
  std::cout << "lambda has been calculated " << std::endl;
  //================================================================
  // Get the mu coefficient

  mu = 0.0;

  // Get -i*omega_c*M[varphi_0]g*

  get_mass(varphi_0_real, g_conj_real, u0, f_real);
  get_mass(varphi_0_real, g_conj_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = -1i * omega_c * (f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  get_mass(varphi_0_imag, g_conj_real, u0, f_imag);
  get_mass(varphi_0_imag, g_conj_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = -1i * omega_c * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  // Get i*omega_c*M[varphi_2]g

  get_mass(varphi_2_real, g_real, u0, f_real);
  get_mass(varphi_2_real, g_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 1i * omega_c * (f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  get_mass(varphi_2_imag, g_real, u0, f_imag);
  get_mass(varphi_2_imag, g_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = 1i * omega_c * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  // Get 0.5*H[u0](g*,varphi_0)

  get_hessian(g_conj_real, varphi_0_real, u0, f_real);
  get_hessian(g_conj_real, varphi_0_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  get_hessian(g_conj_imag, varphi_0_real, u0, f_imag);
  get_hessian(g_conj_imag, varphi_0_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  // Get 0.5*H[u0](g,varphi_2)

  get_hessian(g_real, varphi_2_real, u0, f_real);
  get_hessian(g_real, varphi_2_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  get_hessian(g_imag, varphi_2_real, u0, f_imag);
  get_hessian(g_imag, varphi_2_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  // Get 0.5*H[u0](varphi_0,g*)

  get_hessian(varphi_0_real, g_conj_real, u0, f_real);
  get_hessian(varphi_0_real, g_conj_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  get_hessian(varphi_0_imag, g_conj_real, u0, f_imag);
  get_hessian(varphi_0_imag, g_conj_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  // Get 0.5*H[u0](varphi_2,g)

  get_hessian(varphi_2_real, g_real, u0, f_real);
  get_hessian(varphi_2_real, g_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  get_hessian(varphi_2_imag, g_real, u0, f_imag);
  get_hessian(varphi_2_imag, g_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = 0.5 * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  // Get the K[u0](g,g,g*) + permuatations;
  get_tressian(g_real, g_real, u0, f_real);
  get_tressian(g_real, g_imag, u0, f_imag);

  for (int i = 0; i < N; i++)
  {
    f[i] = (1.0 / 6.0) * (f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  get_tressian(g_imag, g_real, u0, f_imag);
  get_tressian(g_imag, g_imag, u0, f_real);

  for (int i = 0; i < N; i++)
  {
    f[i] = (1.0 / 6.0) * (-1.0 * f_real[i] + f_imag[i]);
  }

  // Update mu
  for (int i = 0; i < N; i++)
  {
    mu += f[i] * g_adjoint[i];
  }

  std::cout << " " << std::endl;
  std::cout << "mu has been calculated " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "===========================" << std::endl;
  std::cout << "ALL COEFFICIENTS CALCULATED" << std::endl;
  std::cout << "===========================" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "nu is " << nu << std::endl;
  std::cout << "lambda is " << lambda << std::endl;
  std::cout << "mu is " << mu << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Scaled by nu " << std::endl;
  std::cout << "lambda_hat is " << lambda / nu << std::endl;
  std::cout << "mu_hat is " << mu / nu << std::endl;
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::fill_in_steady_fourth_order_solution(
  Vector<complex<double>>& varphi_4,
  Vector<complex<double>>& varphi_3,
  int delta)
{
  int N = this->ndof();
  CRDoubleMatrix jacobian(this->dof_distribution_pt());
  CRDoubleMatrix jacobian_0(this->dof_distribution_pt());
  CRDoubleMatrix jacobian_Q(this->dof_distribution_pt());

  CRDoubleMatrix mass(this->dof_distribution_pt());
  DoubleVector residuals(this->dof_distribution_pt());
  DoubleVector g(this->dof_distribution_pt());
  DoubleVector f(this->dof_distribution_pt());
  DoubleVector f_temp(this->dof_distribution_pt());
  Vector<double> u0(N);

  for (int i = 0; i < N; i++)
  {
    g[i] = real(varphi_3[i]);
  }

  for (int i = 0; i < N; i++)
  {
    u0[i] = this->dof(i);
  }

  std::cout << " " << std::endl;
  std::cout << "Creating the rhs for varphi_4..." << std::endl;

  // Here we have the residuals and jacobian for the base problem
  Problem::get_eigenproblem_matrices(mass, jacobian_0);
  Problem::get_jacobian(residuals, jacobian); // This one for the solver

  // Now set Q_inv = 0
  //  set_Q_zero(0);
  Problem::get_jacobian(residuals, jacobian_Q);

  mult_constant(jacobian_Q, -1);
  jacobian_0.add(jacobian_Q, jacobian_0);
  jacobian_0.multiply(g, f_temp);

  //  set_Q_zero(1);

  for (int i = 0; i < N; i++)
  {
    f[i] = -1.0 * delta * f_temp[i];
  }

  for (int i = 0; i < N; i++)
  {
    f_temp[i] = 0.0;
  }

  // Now get hessian term
  get_hessian(g, g, u0, f_temp);

  for (int i = 0; i < N; i++)
  {
    f[i] = (f[i] - 0.5 * f_temp[i]);
  }

  for (int i = 0; i < N; i++)
  {
    f_temp[i] = 0.0;
  }

  std::cout << " " << std::endl;
  std::cout << "Solving for varphi_4..." << std::endl;

  jacobian.solve(f, f_temp);

  ofstream outdata;
  char filename[100];
  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "varphi_4.dat").c_str());
  outdata.open(filename);
  for (int i = 0; i < N; i++)
  {
    varphi_4[i] = f_temp[i] + 1i * 0.0;
    outdata << f_temp[i] << " " << std::endl;
  }

  outdata.close();
  std::cout << " " << std::endl;
  std::cout << "varphi_4 has been calculated." << std::endl;
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::fill_in_second_order_solution(
  Vector<complex<double>>& g,
  double Q_c,
  double omega_c,
  int mode,
  Vector<complex<double>>& varphi_0,
  Vector<complex<double>>& varphi_2,
  Vector<complex<double>>& varphi_3)
{
  int N = this->ndof();
  CRDoubleMatrix jacobian(this->dof_distribution_pt());
  CRDoubleMatrix mass(this->dof_distribution_pt());
  DoubleVector residuals(this->dof_distribution_pt());
  DoubleVector f_real(this->dof_distribution_pt());
  DoubleVector f_imag(this->dof_distribution_pt());
  DoubleVector f_real_temp(this->dof_distribution_pt());
  DoubleVector f_imag_temp(this->dof_distribution_pt());
  DoubleVector g_real(this->dof_distribution_pt());
  DoubleVector g_imag(this->dof_distribution_pt());
  CRDoubleMatrix l_hat_real(this->dof_distribution_pt());
  CRDoubleMatrix l_hat_imag(this->dof_distribution_pt());
  CRDoubleMatrix l_hat_real_temp(this->dof_distribution_pt());
  CRDoubleMatrix l_hat_imag_temp(this->dof_distribution_pt());

  ofstream outdata;
  char filename[100];

  vector<double> u0(N);

  std::cout << "========================================" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Finding the second order solution" << std::endl;

  for (int i = 0; i < N; i++)
  {
    u0[i] = this->dof(i);
  }

  Problem::get_eigenproblem_matrices(mass, jacobian);
  Problem::get_jacobian(residuals, jacobian);

  // Needed for l_hat operator
  Problem::get_eigenproblem_matrices(l_hat_imag, l_hat_real);

  for (int i = 0; i < N; i++)
  {
    g_real[i] = real(g[i]);
    g_imag[i] = imag(g[i]);
  }

  std::cout << " " << std::endl;
  std::cout << "Creating rhs for varphi_0...";

  // Construct M[g_r]g_r + i*M[g_r]g_i;

  get_mass(g_real, g_real, u0, f_real_temp);
  get_mass(g_real, g_imag, u0, f_imag_temp);

  for (int i = 0; i < N; i++)
  {
    f_real[i] += f_real_temp[i];
    f_imag[i] += f_imag_temp[i];
  }

  // Construct -M[g_i]g_i + i*M[g_i]g_r;

  get_mass(g_imag, g_imag, u0, f_real_temp);
  get_mass(g_imag, g_real, u0, f_imag_temp);

  for (int i = 0; i < N; i++)
  {
    f_real[i] -= f_real_temp[i];
    f_imag[i] += f_imag_temp[i];
  }

  // Multiply by -i*omega_c;

  for (int i = 0; i < N; i++)
  {
    f_real[i] = omega_c * f_real[i];
    f_imag[i] = f_real[i];
  }

  // Get Hessian 0.5*H[u0](g,g) = 0.5*H[u0](g_r,g_r) - 0.5*H[u0](g_i,g_i) +
  // 0.5*i*H[u0](g_r,g_i) + 0.5*i*H[u0](g_i,g_r);

  get_hessian(g_real, g_real, u0, f_real_temp);
  get_hessian(g_real, g_imag, u0, f_imag_temp);

  for (int i = 0; i < N; i++)
  {
    f_real[i] += -0.5 * f_real_temp[i];
    f_imag[i] += -0.5 * f_imag_temp[i];
  }

  get_hessian(g_imag, g_real, u0, f_real_temp);
  get_hessian(g_imag, g_real, u0, f_imag_temp);

  for (int i = 0; i < N; i++)
  {
    f_real[i] += 0.5 * f_real_temp[i];
    f_imag[i] += -0.5 * f_imag_temp[i];
  }

  // Create the complex RHS
  Vector<complex<double>> rhs(N);
  Vector<complex<double>> result(N);

  for (int i = 0; i < N; i++)
  {
    rhs[i] = f_real[i] + 1i * f_imag[i];
  }

  std::cout << "..........................done" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Creating l_hat for varphi_0...";

  // Create the complex mass matrix and multiply by mult
  complex<double> mult = complex<double>(2) * 1i * omega_c;

  CRComplexMatrix l_hat;

  // Make 2*i*omega_c
  create_complex_mass_matrix(l_hat, mult);

  // Make 2*i*omega_c*M + J
  l_hat.add(l_hat_real, l_hat);

  std::cout << "............done" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Solving for varphi_0...";

  // Solve the equation
  l_hat.solve(rhs, result);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i];
  }

  for (int i = 0; i < N; i++)
  {
    varphi_0[i] = result[i];
  }

  // Reset stuff
  for (int i = 0; i < N; i++)
  {
    result[i] = 0.0;
  }

  std::cout << "........................done" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "varphi_0 has been calculated" << std::endl;

  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "varphi_0.dat").c_str());

  outdata.open(filename);

  for (int i = 0; i < N; i++)
  {
    outdata << varphi_0[i] << std::endl;
  }

  outdata.close();

  // Reset values
  for (int i = 0; i < N; i++)
  {
    f_real[i] = 0.0;
    f_imag[i] = 0.0;
    f_real_temp[i] = 0.0;
    f_imag_temp[i] = 0.0;
  }

  std::cout << " " << std::endl;
  std::cout << "Creating rhs for varphi_2...";

  // Get M[g](g*) + M[g*](g) =2.0*M[g_real](g_real) + 2.0*M[g_imag](g_imag)
  get_mass(g_real, g_real, u0, f_real_temp);
  get_mass(g_imag, g_imag, u0, f_imag_temp);

  for (int i = 0; i < N; i++)
  {
    f_real[i] = -2.0 * (f_real_temp[i] + f_imag_temp[i]);
    f_imag[i] = 0.0;
  }

  // Multiply by -i*omega_c
  for (int i = 0; i < N; i++)
  {
    f_real[i] = 0.0;
    f_imag[i] = -1.0 * omega_c * f_real[i];
    ;
  }

  // Get Hessian 0.5*H[u0](g,g*) = 0.5*H[u0](g_r,g_r) + 0.5*H[u0](g_i,g_i)
  get_hessian(g_real, g_real, u0, f_real_temp);
  get_hessian(g_imag, g_imag, u0, f_imag_temp);

  for (int i = 0; i < N; i++)
  {
    f_real[i] += -0.5 * (f_real_temp[i] + f_imag_temp[i]);
    f_imag[i] += 0;
  }

  // Create complex rhs
  Vector<complex<double>> f(N);
  for (int i = 0; i < N; i++)
  {
    f[i] = f_real[i] + 1i * f_imag[i];
  }

  std::cout << "...................done" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Solving for varphi_2...";

  // Solve the equation (real and imag parts separately as jacobian is purely
  // real)
  jacobian.solve(f_real, f_real);
  jacobian.solve(f_imag, f_imag);

  for (int i = 0; i < N; i++)
  {
    varphi_2[i] = f_real[i] + 1i * f_imag[i];
  }

  // Reset values
  for (int i = 0; i < N; i++)
  {
    f[i] = 0.0;
    f_real[i] = 0.0;
    f_real_temp[i] = 0.0;
    f_imag[i] = 0.0;
    f_imag_temp[i] = 0.0;
    result[i] = 0.0;
  }

  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "varphi_2.dat").c_str());

  outdata.open(filename);

  for (int i = 0; i < N; i++)
  {
    outdata << real(varphi_2[i]) << " " << imag(varphi_2[i]) << std::endl;
  }

  outdata.close();

  std::cout << "........................done" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "varphi_2 has been calculated" << std::endl;

  // Ideally we set Q_inv = 0.0; instead we say that Q is huge
  set_Q(1e30);

  std::cout << 1.0 / get_Q() << std::endl;

  std::cout << " " << std::endl;
  std::cout << "Creating rhs for varphi_3...";

  // get delta*Q_c*F(u0) term;
  Problem::get_jacobian(f_real_temp, jacobian);

  for (int i = 0; i < N; i++)
  {
    f_real_temp[i] = residuals[i] - f_real_temp[i];
  }

  // Reset Q
  set_Q(Q_c);
  Problem::get_jacobian(residuals, jacobian);

  std::cout << "...................done" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Solving for varphi_3...";

  jacobian.solve(f_real_temp, residuals);

  for (int i = 0; i < N; i++)
  {
    varphi_3[i] = residuals[i] + 1i * 0.0;
    f_real_temp[i] = 0.0;
  }

  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "varphi_3.dat").c_str());

  outdata.open(filename);

  for (int i = 0; i < N; i++)
  {
    outdata << real(varphi_3[i]) << " " << imag(varphi_3[i]) << std::endl;
  }

  outdata.close();

  std::cout << "............done" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "varphi_3 has been calculated" << std::endl;

  std::cout << " " << std::endl;
  std::cout << "=============================================" << std::endl;
  std::cout << "THE SECOND ORDER SOLUTION HAS BEEN CALCULATED" << std::endl;
  std::cout << "=============================================" << std::endl;
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::get_mass(DoubleVector& g1,
                                               DoubleVector& g2,
                                               vector<double>& u0,
                                               DoubleVector& f)
{
  int N = this->ndof();
  CRDoubleMatrix mass_plus(this->dof_distribution_pt());
  CRDoubleMatrix mass_neg(this->dof_distribution_pt());
  CRDoubleMatrix jacobian(this->dof_distribution_pt());
  DoubleVector residuals(this->dof_distribution_pt());
  double epsilon = 1e-8;

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i] + epsilon * g1[i];
  }

  Problem::get_eigenproblem_matrices(mass_plus, jacobian);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i] - epsilon * g1[i];
  }

  Problem::get_eigenproblem_matrices(mass_neg, jacobian);

  mult_constant(mass_neg, -1.0);

  mass_plus.add(mass_neg, jacobian);

  mult_constant(jacobian, -0.5 / epsilon);

  jacobian.multiply(g2, f);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i];
  }
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::get_hessian(DoubleVector& g1,
                                                  DoubleVector& g2,
                                                  vector<double>& u0,
                                                  DoubleVector& f)
{
  int N = this->ndof();
  CRDoubleMatrix jacobian_plus(this->dof_distribution_pt());
  CRDoubleMatrix jacobian_neg(this->dof_distribution_pt());
  DoubleVector residuals(this->dof_distribution_pt());

  double epsilon = 1e-7;

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i] + epsilon * g1[i];
  }

  Problem::get_jacobian(residuals, jacobian_plus);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i] - epsilon * g1[i];
  }

  Problem::get_jacobian(residuals, jacobian_neg);

  // Calculate (j_+ - j_-)/(2*eps);
  mult_constant(jacobian_neg, -1.0);

  jacobian_plus.add(jacobian_neg, jacobian_plus);

  double mult = 1.0 / (2 * epsilon);

  mult_constant(jacobian_plus, mult);

  // Calculate the matrix/vector product
  jacobian_plus.multiply(g2, f);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i];
  }
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::get_tressian(DoubleVector& g1,
                                                   DoubleVector& g2,
                                                   vector<double>& u0,
                                                   DoubleVector& f)
{
  int N = this->ndof();
  CRDoubleMatrix jacobian_plus(this->dof_distribution_pt());
  CRDoubleMatrix jacobian_neg(this->dof_distribution_pt());
  CRDoubleMatrix jacobian_0(this->dof_distribution_pt());
  DoubleVector residuals(this->dof_distribution_pt());

  double epsilon = 1e-4;

  Problem::get_jacobian(residuals, jacobian_0);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i] + epsilon * g1[i];
  }

  Problem::get_jacobian(residuals, jacobian_plus);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i] - epsilon * g1[i];
  }

  Problem::get_jacobian(residuals, jacobian_neg);

  // Calculate (j_+ -2*j_0  + j_-)/(eps*eps);
  mult_constant(jacobian_0, -2.0);

  jacobian_plus.add(jacobian_neg, jacobian_plus);
  jacobian_plus.add(jacobian_0, jacobian_plus);

  double mult = 1.0 / (epsilon * epsilon);

  mult_constant(jacobian_plus, mult);

  // Calculate the matrix/vector product
  jacobian_plus.multiply(g2, f);

  for (int i = 0; i < N; i++)
  {
    this->dof(i) = u0[i];
  }
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::mult_constant(CRComplexMatrix& A,
                                                    complex<double> mult)
{
  CRComplexMatrix identity;
  unsigned nrow_local = identity.nrow();
  unsigned first_row = 0;

  // Diagonal Entries
  Vector<std::complex<double>> values(nrow_local, mult);

  // Column indices linked to values
  Vector<int> column_indices(nrow_local);
  Vector<int> row_start(nrow_local + 1);

  for (unsigned i = 0; i < nrow_local; ++i)
  {
    column_indices[i] = first_row + i;
    row_start[i] = i;
  }

  row_start[nrow_local] = nrow_local;

  // Build

  unsigned ncol = this->ndof();

  identity.build(values, column_indices, row_start, ncol, ncol);

  std::cout << "I am here" << std::endl;

  int nnz = A.nnz();

  int* col = A.column_index();
  int* row = A.row_start();

  ofstream outdata;
  char filename[100];

  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "Rgrifter1.dat").c_str());
  outdata.open(filename);

  A.sparse_indexed_output_helper(outdata);

  outdata.close();

  int N = this->ndof();

  for (unsigned long i = 0; i < 40; i++)
  {
    for (long j = row[i]; j < row[i + 1]; j++)
    {
      std::cout << i << " " << col[j] << " " << A(i, col[j]) << std::endl;
    }
  }

  std::cout << "I am here2" << std::endl;

  std::cout << row[0] << std::endl;
  std::cout << row[1] << std::endl;
  std::cout << row[2] << std::endl;
  std::cout << row[nnz + 1] << std::endl;

  A.multiply(identity, A);
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::mult_constant(CRDoubleMatrix& A,
                                                    double mult)
{
  CRDoubleMatrix identity(this->dof_distribution_pt());
  unsigned nrow_local = identity.distribution_pt()->nrow_local();
  unsigned first_row = identity.distribution_pt()->first_row();

  // Diagonal Entries
  Vector<double> values(nrow_local, mult);

  // Column indices linked to values
  Vector<int> column_indices(nrow_local);
  Vector<int> row_start(nrow_local + 1);

  for (unsigned i = 0; i < nrow_local; ++i)
  {
    column_indices[i] = first_row + i;
    row_start[i] = i;
  }

  row_start[nrow_local] = nrow_local;

  // Build

  unsigned ncol = this->ndof();

  identity.build(ncol, values, column_indices, row_start);

  // Multiply

  A.multiply(identity, A);
}


template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::test_identity_matrix()
{
  int N = this->ndof();
  CRDoubleMatrix jaco(this->dof_distribution_pt());
  CRDoubleMatrix mass(this->dof_distribution_pt());
  CRDoubleMatrix jaco_transpose(this->dof_distribution_pt());
  jaco_transpose.build(this->dof_distribution_pt());

  DoubleVector residuals(this->dof_distribution_pt());
  DoubleVector f(this->dof_distribution_pt());

  Problem::get_jacobian(residuals, jaco);

  CRDoubleMatrix result1;
  result1.build(this->dof_distribution_pt());

  CRDoubleMatrix result2;
  result2.build(this->dof_distribution_pt());

  ////////////////////////////////////////////////////////
  /// MAKE A CONSTANT DIAGONAL MATRIX//////////////////////

  CRDoubleMatrix ident;
  ident.build(this->dof_distribution_pt());

  unsigned nrow_local = ident.distribution_pt()->nrow_local();
  unsigned first_row = ident.distribution_pt()->first_row();

  // Diagonal Entries
  Vector<double> values(nrow_local, 2.0);

  // Column indices linked to values
  Vector<int> column_indices(nrow_local);
  Vector<int> row_start(nrow_local + 1);

  for (unsigned i = 0; i < nrow_local; ++i)
  {
    column_indices[i] = first_row + i;
    row_start[i] = i;
  }

  row_start[nrow_local] = nrow_local;

  // Build

  unsigned ncol = this->ndof();

  ident.build(ncol, values, column_indices, row_start);

  std::cout << ident(0, 0) << " ";
  /////////////////////////////////////////////////////////

  CRDoubleMatrix result3;
  ident.multiply(jaco, result3);

  jaco.multiply(ident, result1);
  jaco.add(jaco, result2);

  std::cout << "Number of unknowns " << N << std::endl;
  std::cout << "Number of non-zeros entries of J " << jaco.nnz() << std::endl;
  std::cout << "Number of non-zeros entries of I " << ident.nnz() << std::endl;
  std::cout << "Number of non-zeros entries of J*I " << result1.nnz()
            << std::endl;
  std::cout << "Number of non-zeros entries of J+J " << result2.nnz()
            << std::endl;
  std::cout << "Number of non-zeroes of I*J " << result3.nnz() << std::endl;
  std::cout << " " << std::endl;
  std::cout << "J(0,107) entry is " << jaco(0, 107) << std::endl;
  std::cout << "J(107,0) entry is " << jaco(107, 0) << std::endl;

  std::cout << "J*I(0,107) entry is " << result1(0, 107) << std::endl;
  std::cout << "[J+J](0,107) entry is " << result2(0, 107) << std::endl;

  // J_T.build(this->dof_distribution_pt);
  CRDoubleMatrix* J_T_pt = &jaco_transpose;

  jaco.get_matrix_transpose(J_T_pt);

  std::cout << "J_T(0,107) entry is " << jaco_transpose(0, 107) << std::endl;
  std::cout << "J_T(107,0) entry is " << jaco_transpose(107, 0) << std::endl;
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_complex_mass_matrix(
  CRComplexMatrix& mass_complex, complex<double> mult)
{
  // Create the space for the mass and jacobian matrices
  CRDoubleMatrix mass(this->dof_distribution_pt());
  CRDoubleMatrix jacobian(this->dof_distribution_pt());

  // Build the matrices
  Problem::get_eigenproblem_matrices(mass, jacobian);

  // The number of non-zero entries for the mass matrix
  unsigned long nnz = mass.nnz();

  // Vectors for the CRDouble mass matrix
  double* values = mass.value();
  unsigned N = mass.nrow_local();

  // Output the CRDouble mass
  mass.sparse_indexed_output_with_offset(
    Problem_Parameter::Doc_info.directory() + "sparse_mass.dat");

  // Output the CRDouble jacobian
  jacobian.sparse_indexed_output_with_offset(
    Problem_Parameter::Doc_info.directory() + "sparse_jacobian.dat");

  // Create vector of values and the pointer
  Vector<complex<double>> value(nnz);

  for (unsigned i = 0; i < nnz; i++)
  {
    value[i] = mult * values[i];
  }

  // Build the CRComplex column and row vectors directly from the known
  // distribution of the CRDoubleMatrix
  Vector<int> col(nnz);
  Vector<int> row(N + 1);

  for (unsigned i = 0; i < N; i++)
  {
    row[i] = mass.row_start()[i];
  }

  for (unsigned i = 0; i < nnz; i++)
  {
    col[i] = mass.column_index()[i];
  }

  // Build the CRComplex mass matrix from vectors
  mass_complex.build(value, col, row, N, N);

  // Output the CRComplexMatrix
  ofstream outdata;
  char filename[100];
  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "sparse_complex_mass.dat")
            .c_str());
  outdata.open(filename);
  for (long int i = 0; i < N; i++)
  {
    for (int j = row[i]; j < row[i + 1]; j++)
    {
      outdata << i << " " << col[j] << " " << real(value[i]) << " "
              << imag(value[j]) << std::endl;
    }
  }
  outdata.close();
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::create_complex_jacobian_matrix(
  CRComplexMatrix& jacobian_complex, complex<double> mult)
{
  // Create the space for the mass and jacobian matrices
  CRDoubleMatrix mass(this->dof_distribution_pt());
  CRDoubleMatrix jacobian(this->dof_distribution_pt());

  // Build the matrices
  Problem::get_eigenproblem_matrices(mass, jacobian);

  // The number of non-zero entries for the mass matrix
  unsigned long nnz = jacobian.nnz();

  // Vectors for the CRDouble mass matrix
  double* values = jacobian.value();
  unsigned N = jacobian.nrow_local();

  // Output the CRDouble mass
  mass.sparse_indexed_output_with_offset(
    Problem_Parameter::Doc_info.directory() + "sparse_mass.dat");

  // Output the CRDouble jacobian
  jacobian.sparse_indexed_output_with_offset(
    Problem_Parameter::Doc_info.directory() + "sparse_jacobian.dat");

  // Create vector of values and the pointer
  Vector<complex<double>> value(nnz);

  for (unsigned i = 0; i < nnz; i++)
  {
    value[i] = mult * values[i];
  }

  // Build the CRComplex column and row vectors directly from the known
  // distribution of the CRDoubleMatrix
  Vector<int> col(nnz);
  Vector<int> row(N + 1);

  for (unsigned i = 0; i < N; i++)
  {
    row[i] = jacobian.row_start()[i];
  }

  for (unsigned i = 0; i < nnz; i++)
  {
    col[i] = jacobian.column_index()[i];
  }

  // Build the CRComplex mass matrix from vectors
  jacobian_complex.build(value, col, row, N, N);

  // Output the CRComplexMatrix
  ofstream outdata;
  char filename[100];
  sprintf(
    filename,
    (Problem_Parameter::Doc_info.directory() + "sparse_complex_jacobian.dat")
      .c_str());
  outdata.open(filename);
  for (long int i = 0; i < N; i++)
  {
    for (int j = row[i]; j < row[i + 1]; j++)
    {
      outdata << i << " " << col[j] << " " << real(value[i]) << " "
              << imag(value[j]) << std::endl;
    }
  }
  outdata.close();
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::test_complex_solve()
{
  int N = 2;
  int nnz = 4;
  Vector<int> row(N + 1);
  Vector<int> col(nnz);
  Vector<complex<double>> val(nnz);
  Vector<complex<double>> rhs(N);
  Vector<complex<double>> result(N);

  row[0] = 0;
  row[1] = 1;

  col[0] = 0;
  col[1] = 1;
  col[2] = 0;
  col[3] = 1;

  val[0] = 1i;
  val[1] = 2.0;
  val[2] = 4.0;
  val[3] = -1i;

  CRComplexMatrix test;
  test.build(val, col, row, N, N);

  rhs[0] = 1.0 - 2.0 * 1i;
  rhs[1] = -1.0 + 3.0 * 1i;

  cout << "test.solve(rhs, result);" << endl;
  test.solve(rhs, result);

  std::cout << "x is " << result[0] << std::endl;
  std::cout << "y is " << result[1] << std::endl;
}


template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::test_complex_identity_matrix()
{
  // Create the space for the mass and jacobian matrices
  CRDoubleMatrix mass(this->dof_distribution_pt());
  CRDoubleMatrix jacobian(this->dof_distribution_pt());

  // Build the matrices
  Problem::get_eigenproblem_matrices(mass, jacobian);

  // The number of non-zero entries for the mass matrix
  unsigned long nnz = mass.nnz();

  // The multiplier
  complex<double> mult = complex<double>(2) * 1i;

  // Vectors for the CRDouble mass matrix
  double* values = mass.value();
  unsigned N = mass.nrow_local();

  // Output the CRDouble mass
  mass.sparse_indexed_output_with_offset(
    Problem_Parameter::Doc_info.directory() + "sparse_mass.dat");

  // Create vector of values
  Vector<complex<double>> value(nnz);

  for (unsigned i = 0; i < nnz; i++)
  {
    value[i] = mult * values[i];
  }

  // Build the CRComplex column and row vectors directly from the known
  // distribution of the CRDoubleMatrix mass
  Vector<int> col(nnz);
  Vector<int> row(N + 1);

  for (unsigned i = 0; i < N; i++)
  {
    row[i] = mass.row_start()[i];
  }

  for (unsigned i = 0; i < nnz; i++)
  {
    col[i] = mass.column_index()[i];
  }

  // Build the CRComplex mass matrix from vectors
  CRComplexMatrix mass_complex;
  mass_complex.build(value, col, row, N, N);

  // Output the CRComplexMatrix
  ofstream outdata;
  char filename[100];
  sprintf(filename,
          (Problem_Parameter::Doc_info.directory() + "sparse_complex_mass.dat")
            .c_str());
  outdata.open(filename);
  mass_complex.sparse_indexed_output_helper(outdata);
  outdata.close();

  cout << "mass_complex.multiply(mass_complex,mass_complex);" << endl;
  mass_complex.multiply(mass_complex, mass_complex);
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::dump_it(ofstream& dump_file)
{
  // Call generic dump
  dump(dump_file);
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::restart(ifstream& restart_file)
{
  read(restart_file);
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::set_initial_condition()
{
  ifstream* restart_file_pt =
    new ifstream(Problem_Parameter::restart_input_filename, ios_base::in);
  restart(*restart_file_pt);
}

template<class ELEMENT>
void BubbleInChannelProblem<ELEMENT>::set_initial_condition_asymm()
{
  ifstream* restart_file_pt =
    new ifstream(Problem_Parameter::restart_input_filename, ios_base::in);
  restart(*restart_file_pt);
}
#endif
