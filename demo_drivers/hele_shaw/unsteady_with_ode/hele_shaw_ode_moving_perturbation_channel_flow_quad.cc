#include <cmath>
#include <iostream>

#include "generic.h"
#include "meshes.h"
#include "hele_shaw.h"
#include "ode.h"

using namespace oomph;
using namespace std;

namespace problem_parameter {

double *global_time_pt = 0;
double global_reference_frame_speed = 0.0;

void upper_wall_fct(const Vector<double> &x, double &b, double &dbdt) {

  double height = 0.5;
  double rms_width = 0.1;
  double centre_x = 1.0 + (*global_time_pt) * (global_reference_frame_speed);
  double centre_y = 0.5;

  // Transform y such that the domain is between 0 and 1 rather than -1 and 1
  double local_x = x[0];
  double local_y = x[1];

  double f = 0.0;
  f = -(local_x - centre_x) * (local_x - centre_x) /
          (2.0 * rms_width * rms_width) -
      (local_y - centre_y) * (local_y - centre_y) /
          (2.0 * rms_width * rms_width);
  double dfdt = 0.0;
  dfdt = -2.0 * global_reference_frame_speed * (local_x - centre_x) /
         (2.0 * rms_width * rms_width);

  b = 1.0 - height * std::exp(f);

  dbdt = -dfdt * height * std::exp(f);
}

void get_dirichlet_bc(const Vector<double> &x, double &value) { value = 0.0; }

void get_neumann_bc(const Vector<double> &x, double &flux) {
  /// Ahead of the bubble we have the boundary condition dp/dx = -G.
  /// We need to supply the function b^3 n.grad p = -G b^3.
  double b;
  double dbdt;
  upper_wall_fct(x, b, dbdt);
  double G = 1.0;
  double p_x = -1.0 * G;
  flux = p_x * (b * b * b);
}

class DistanceODE : public SolutionFunctorBase {
public:
  DistanceODE() {}

  virtual ~DistanceODE() {}

  Vector<double> operator()(const double &t, const Vector<double> &x) const {
    Vector<double> output(1);
    output[0] = 0.0;
    return output;
  }

  Vector<double> derivative(const double &t, const Vector<double> &x,
                            const Vector<double> &u) const {
    Vector<double> output(1);
    output[0] = global_reference_frame_speed;
    return output;
  }
};

} // namespace problem_parameter

template <class ELEMENT> class HeleShawChannelProblem : public Problem {

public:
  /// Constructor
  HeleShawChannelProblem();

  /// Destructor (empty)
  ~HeleShawChannelProblem() {}

  void set_initial_condition();

  /// Iterate forward in time
  void iterate_timestepper(const double &t_step, const double &t_final,
                           DocInfo &doc_info);

  /// Doc the solution
  void doc_solution(DocInfo &doc_info);

private:
  /// Generate mesh
  void generate_mesh();

  /// Create flux elements
  void create_flux_elements(const unsigned &boundary);

  /// Generate mesh
  void assign_mesh();

  /// Pin dirichlet outlet boundary
  void pin_dirichlet_boundaries();

  /// Upcast elements and finalise setup
  void setup_elements();

  /// Set boundary condition values
  void set_boundary_conditions();

  /// Update the problem specs before solve
  void actions_before_newton_solve();

  /// Update the problem specs before solve (empty)
  void actions_after_newton_solve();

  /// Update the problem specs before timestep
  void actions_before_implicit_timestep();

  /// Save the boundary data to file
  void save_boundaries_to_file(ofstream &output_stream, string filename);

  /// Save the solution data to file
  void save_solution_to_file(ofstream &output_stream, string filename,
                             unsigned n_points);

  /// Save the distance data to file
  void save_distance_to_file(ofstream &output_stream, string filename);

  /// Pointer to the "bulk" mesh
  SimpleRectangularQuadMesh<ELEMENT> *Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh *Surface_mesh_pt;

  /// Pointer to the "distance" mesh
  Mesh *Distance_mesh_pt;

  /// Pointer to the ode element
  ODEElement *ode_element_pt;
};

template <class ELEMENT>
HeleShawChannelProblem<ELEMENT>::HeleShawChannelProblem() {
  cout << "Problem constructor" << endl;

  this->add_time_stepper_pt(new BDF<1>);

  this->generate_mesh();

  this->assign_mesh();

  this->pin_dirichlet_boundaries();

  this->setup_elements();

  // Setup equation numbering scheme
  cout << "Assign equation numbers." << endl;
  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::doc_solution(DocInfo &doc_info) {
  cout << "Doc solution" << endl;

  string data_directory = doc_info.directory();
  ofstream output_stream;

  string filename = data_directory + "/boundaries.dat";
  this->save_boundaries_to_file(output_stream, filename);

  unsigned n_points = 5;
  this->save_solution_to_file(output_stream,
                              data_directory + "/soln" +
                                  to_string(doc_info.number()) + ".dat",
                              n_points);

  this->save_distance_to_file(output_stream, data_directory + "/distance.dat");
}

template <class ELEMENT> void HeleShawChannelProblem<ELEMENT>::generate_mesh() {
  cout << "Generate mesh" << endl;
  unsigned n_x = 32;
  unsigned n_y = 32;
  double l_x = 2.0;
  double l_y = 1.0;

  this->Bulk_mesh_pt = new SimpleRectangularQuadMesh<ELEMENT>(
      n_x, n_y, l_x, l_y, this->time_stepper_pt());

  this->Surface_mesh_pt = new Mesh;

  cout << "Create flux elements" << endl;
  const unsigned flux_boundary = 3;
  this->create_flux_elements(flux_boundary);

  this->Distance_mesh_pt = new Mesh;

  cout << "Creating DistanceODE elements" << endl;
  SolutionFunctorBase *exact_solution_pt = new problem_parameter::DistanceODE;
  ode_element_pt = new ODEElement(this->time_stepper_pt(), exact_solution_pt);

  this->Distance_mesh_pt->add_element_pt(ode_element_pt);
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::create_flux_elements(
    const unsigned &boundary) {

  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(boundary);
  for (unsigned n = 0; n < n_element; n++) {
    ELEMENT *bulk_element_pt = dynamic_cast<ELEMENT *>(
        this->Bulk_mesh_pt->boundary_element_pt(boundary, n));

    int face_index = this->Bulk_mesh_pt->face_index_at_boundary(boundary, n);

    HeleShawFluxElement<ELEMENT> *flux_element_pt =
        new HeleShawFluxElement<ELEMENT>(bulk_element_pt, face_index);

    this->Surface_mesh_pt->add_element_pt(flux_element_pt);
  }
}

template <class ELEMENT> void HeleShawChannelProblem<ELEMENT>::assign_mesh() {
  cout << "Assign mesh" << endl;
  this->add_sub_mesh(this->Bulk_mesh_pt);
  this->add_sub_mesh(this->Surface_mesh_pt);
  this->add_sub_mesh(this->Distance_mesh_pt);

  this->build_global_mesh();
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::pin_dirichlet_boundaries() {
  cout << "Pin Dirichlet boundaries" << endl;
  unsigned n_boundary = this->Bulk_mesh_pt->nboundary();
  bool pin_boundary[n_boundary] = {false};
  pin_boundary[1] = true;
  for (unsigned b = 0; b < n_boundary; b++) {
    if (pin_boundary[b]) {
      cout << "Pinning boundary: " << b << endl;
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_node; n++) {
        this->Bulk_mesh_pt->boundary_node_pt(b, n)->pin(0);
      }
    }
  }
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::setup_elements() {
  cout << "Setup elements" << endl;

  // Find number of elements in mesh
  unsigned n_element = this->Bulk_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++) {
    // Upcast from GeneralElement to the present element
    ELEMENT *el_pt = dynamic_cast<ELEMENT *>(this->Bulk_mesh_pt->element_pt(i));

    el_pt->upper_wall_fct_pt() = problem_parameter::upper_wall_fct;
  }

  // Find number of elements in mesh
  n_element = this->Surface_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++) {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElement<ELEMENT> *el_pt =
        dynamic_cast<HeleShawFluxElement<ELEMENT> *>(
            this->Surface_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = &problem_parameter::get_neumann_bc;
  }
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::set_boundary_conditions() {

  unsigned n_boundary = this->Bulk_mesh_pt->nboundary();
  bool set_boundary[n_boundary] = {false};
  set_boundary[1] = true;
  for (unsigned b = 0; b < n_boundary; b++) {
    if (set_boundary[b]) {
      cout << "Setting boundary: " << b << endl;
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_node; n++) {
        Node *node_pt = this->Bulk_mesh_pt->boundary_node_pt(b, n);
        double value = 0.0;
        Vector<double> x(2);
        x[0] = node_pt->x(0);
        x[1] = node_pt->x(1);
        problem_parameter::get_dirichlet_bc(x, value);
        node_pt->set_value(0, value);
      }
    }
  }
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::set_initial_condition() {
  unsigned element_index = 0;
  unsigned data_index = 0;
  unsigned time_level = 0;
  unsigned value_index = 0;
  double x0 = 0.0;
  this->Distance_mesh_pt->element_pt(element_index)
      ->internal_data_pt(data_index)
      ->set_value(time_level, value_index, x0);
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::actions_before_newton_solve() {
  this->set_boundary_conditions();

  unsigned element_index = 0;
  unsigned data_index = 0;
  unsigned value_index = 0;
  this->Distance_mesh_pt->element_pt(element_index)
      ->internal_data_pt(data_index)
      ->pin(value_index);

  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;

  cout << this->self_test() << endl;
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::actions_after_newton_solve() {
  unsigned element_index = 0;
  unsigned data_index = 0;
  unsigned value_index = 0;
  this->Distance_mesh_pt->element_pt(element_index)
      ->internal_data_pt(data_index)
      ->unpin(value_index);

  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::actions_before_implicit_timestep() {}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::iterate_timestepper(const double &t_step,
                                                          const double &t_final,
                                                          DocInfo &doc_info) {
  unsigned n_timestep = ceil(t_final / t_step);

  for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++) {
    cout << "t: " << this->time_pt()->time() << endl;

    this->unsteady_newton_solve(t_step);

    doc_info.number()++;
    this->doc_solution(doc_info);
  }
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_boundaries_to_file(
    ofstream &output_stream, string filename) {
  output_stream.open(filename.c_str());
  this->Bulk_mesh_pt->output_boundaries(output_stream);
  output_stream.close();
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_solution_to_file(
    ofstream &output_stream, string filename, unsigned n_points) {
  output_stream.open(filename);
  this->Bulk_mesh_pt->output(output_stream, n_points);
  output_stream.close();
}

template <class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_distance_to_file(
    ofstream &output_stream, string filename) {
  output_stream.open(filename);
  unsigned element_index = 0;
  unsigned data_index = 0;
  double t = this->Distance_mesh_pt->element_pt(element_index)
                 ->internal_data_pt(data_index)
                 ->time_stepper_pt()
                 ->time();
  unsigned value_index = 0;
  double x = this->Distance_mesh_pt->element_pt(element_index)
                 ->internal_data_pt(data_index)
                 ->value(value_index);
  output_stream << "t: " << t << ", x: " << x << endl;
  output_stream.close();
}

int main(int argc, char *argv[]) {
  cout << "Hele-Shaw channel flow" << endl;

  /// Store command line arguments
  CommandLineArgs::setup(argc, argv);

  /// Label for output
  DocInfo doc_info;

  /// Output directory
  doc_info.set_directory("data/");

  HeleShawChannelProblem<QHeleShawElement<3>> problem;

  cout << "\n\n\nProblem self-test ";
  if (problem.self_test() == 0) {
    cout << "passed: Problem can be solved." << std::endl;
  } else {
    throw OomphLibError("Self test failed", OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  problem_parameter::global_time_pt = &problem.time_pt()->time();
  problem_parameter::global_reference_frame_speed = 1.0;

  // Solve for initial conditions
  cout << "Newton solve" << endl;
  problem.newton_solve();
  cout << "t: " << problem.time_pt()->time() << endl;

  double dt = 0.1;
  double tF = 0.5;
  cout << "Iterate timestepper solve" << endl;
  problem.iterate_timestepper(dt, tF, doc_info);

  cout << "End of Hele-Shaw channel flow" << endl;

  return 0;
}
