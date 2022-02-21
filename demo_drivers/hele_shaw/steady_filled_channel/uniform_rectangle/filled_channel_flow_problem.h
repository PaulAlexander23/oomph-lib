#ifndef OOMPH_FILLED_CHANNEL_FLOW_PROBLEM_HEADER
#define OOMPH_FILLED_CHANNEL_FLOW_PROBLEM_HEADER

#include "generic.h"
//#include "meshes.h"
#include "hele_shaw.h"
#include "info_element.h"
#include "problem_parameter.h"

template<class ELEMENT>
class HeleShawChannelProblem : public Problem
{
public:
  /// Constructor
  HeleShawChannelProblem() {}

  /// Destructor (empty)
  ~HeleShawChannelProblem() {}


  void setup();

  /// Doc the solution
  void doc_solution(DocInfo& doc_info);

  /// Generate bulk mesh
  virtual void generate_bulk_mesh() = 0;

protected:
  /// Pointer to the "bulk" mesh
  Mesh* Bulk_mesh_pt;

private:
  /// Generate mesh
  void generate_mesh();

  /// Generate info mesh
  void generate_info_mesh();

  /// Generate flux elements
  void generate_inlet_flux_mesh();
  void generate_outlet_flux_mesh();

  /// Pin dirichlet outlet boundary
  void pin_data();

  /// Upcast elements and finalise setup
  void upcast_and_finalise_elements();

  /// Set boundary condition values
  void set_boundary_conditions();

  /// Save the boundary data to file
  void save_boundaries_to_file(ofstream& output_stream, string filename);

  /// Save the solution data to file
  void save_solution_to_file(ofstream& output_stream,
                             string filename,
                             unsigned n_points);

  /// Save the integral to file
  void save_integral_to_file(ofstream& output_stream, string filename);

  /// Pointer to the "flux" mesh
  Mesh* Inlet_flux_mesh_pt;
  Mesh* Outlet_flux_mesh_pt;

  /// Pointer to the info mesh
  Mesh* Info_mesh_pt;

  /// Pointer to the integral data
  Data* Inlet_integral_data_pt;
  Data* Outlet_integral_data_pt;
};

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::setup()
{
  cout << "generate_mesh" << endl;
  this->generate_mesh();

  cout << "upcast_and_finalise_elements" << endl;
  this->upcast_and_finalise_elements();

  cout << "pin_data" << endl;
  this->pin_data();

  cout << "set_boundary_conditions" << endl;
  this->set_boundary_conditions();

  // Setup equation numbering scheme
  cout << "Number of equations: " << this->assign_eqn_numbers() << endl;
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  string data_directory = doc_info.directory();
  ofstream output_stream;

  string filename = data_directory + "/boundaries.dat";
  this->save_boundaries_to_file(output_stream, filename);

  unsigned n_points = 5;
  this->save_solution_to_file(
    output_stream, data_directory + "/soln.dat", n_points);

  filename = data_directory + "/integral.dat";
  this->save_integral_to_file(output_stream, filename);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_mesh()
{
  Inlet_integral_data_pt = new Data(1);
  Outlet_integral_data_pt = new Data(1);

  generate_bulk_mesh();
  cout << "finished generate_bulk_mesh" << endl;
  this->generate_info_mesh();
  this->generate_inlet_flux_mesh();
  this->generate_outlet_flux_mesh();

  this->add_sub_mesh(this->Bulk_mesh_pt);
  this->add_sub_mesh(this->Inlet_flux_mesh_pt);
  this->add_sub_mesh(this->Outlet_flux_mesh_pt);
  this->add_sub_mesh(this->Info_mesh_pt);

  this->build_global_mesh();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_info_mesh()
{
  this->Info_mesh_pt = new Mesh;

  this->Info_mesh_pt->add_element_pt(new InfoElement(Inlet_integral_data_pt));
  this->Info_mesh_pt->add_element_pt(new InfoElement(Outlet_integral_data_pt));
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_inlet_flux_mesh()
{
  const unsigned boundary = 3;

  this->Inlet_flux_mesh_pt = new Mesh;

  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(boundary);
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(boundary, n));

    int face_index = this->Bulk_mesh_pt->face_index_at_boundary(boundary, n);

    HeleShawFluxElement<ELEMENT>* flux_element_pt =
      new HeleShawFluxElement<ELEMENT>(bulk_element_pt, face_index);

    flux_element_pt->add_b3_data_pt(
      this->Info_mesh_pt->element_pt(0)->internal_data_pt(0));

    this->Inlet_flux_mesh_pt->add_element_pt(flux_element_pt);
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::generate_outlet_flux_mesh()
{
  const unsigned boundary = 1;

  this->Outlet_flux_mesh_pt = new Mesh;

  unsigned n_element = this->Bulk_mesh_pt->nboundary_element(boundary);
  for (unsigned n = 0; n < n_element; n++)
  {
    ELEMENT* bulk_element_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(boundary, n));

    int face_index = this->Bulk_mesh_pt->face_index_at_boundary(boundary, n);

    HeleShawFluxElement<ELEMENT>* flux_element_pt =
      new HeleShawFluxElement<ELEMENT>(bulk_element_pt, face_index);

    flux_element_pt->add_b3_data_pt(
      this->Info_mesh_pt->element_pt(1)->internal_data_pt(0));

    this->Outlet_flux_mesh_pt->add_element_pt(flux_element_pt);
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::pin_data()
{
  unsigned n_boundary = this->Bulk_mesh_pt->nboundary();
  bool pin_boundary[n_boundary] = {false};
  pin_boundary[1] = true;
  for (unsigned b = 0; b < n_boundary; b++)
  {
    if (pin_boundary[b])
    {
      unsigned n_node = this->Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n = 0; n < n_node; n++)
      {
        this->Bulk_mesh_pt->boundary_node_pt(b, n)->pin(0);
      }
    }
  }

  unsigned index = 0;
  this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->unpin(index);
  this->Info_mesh_pt->element_pt(1)->internal_data_pt(0)->unpin(index);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::upcast_and_finalise_elements()
{
  /// Bulk mesh
  // Find number of elements in mesh
  unsigned n_element = this->Bulk_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->Bulk_mesh_pt->element_pt(i));

    el_pt->upper_wall_fct_pt() = problem_parameter::upper_wall_fct;
  }

  /// Info mesh
  unsigned index = 0;
  problem_parameter::inlet_area_pt =
    this->Info_mesh_pt->element_pt(0)->internal_data_pt(0)->value_pt(index);
  problem_parameter::outlet_area_pt =
    this->Info_mesh_pt->element_pt(1)->internal_data_pt(0)->value_pt(index);

  /// Inlet mesh
  // Find number of elements in mesh
  n_element = this->Inlet_flux_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElement<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElement<ELEMENT>*>(
        this->Inlet_flux_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = &problem_parameter::get_inlet_flux_bc;
  }

  /// Outlet mesh
  // Find number of elements in mesh
  n_element = this->Outlet_flux_mesh_pt->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for (unsigned i = 0; i < n_element; i++)
  {
    // Upcast from GeneralElement to the present element
    HeleShawFluxElement<ELEMENT>* el_pt =
      dynamic_cast<HeleShawFluxElement<ELEMENT>*>(
        this->Outlet_flux_mesh_pt->element_pt(i));

    // Set the Neumann function pointer
    el_pt->flux_fct_pt() = &problem_parameter::get_outlet_flux_bc;
  }
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::set_boundary_conditions()
{
  /// Pin a single node in the bulk mesh
  /// Pinning the zeroth node on the zeroth boundary
  const unsigned b = 0;
  const unsigned n = 0;
  Node* node_pt = this->Bulk_mesh_pt->boundary_node_pt(b, n);
  double value = 0.0;
  Vector<double> x(2);
  x[0] = node_pt->x(0);
  x[1] = node_pt->x(1);
  problem_parameter::get_dirichlet_bc(x, value);
  node_pt->set_value(0, value);

  cout << "Integral info mesh" << endl;
  /// Integral info mesh
  Inlet_integral_data_pt->set_value(0, 1.0);
  Outlet_integral_data_pt->set_value(0, 1.0);
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_boundaries_to_file(
  ofstream& output_stream, string filename)
{
  output_stream.open(filename.c_str());
  this->Bulk_mesh_pt->output_boundaries(output_stream);
  output_stream.close();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_solution_to_file(
  ofstream& output_stream, string filename, unsigned n_points)
{
  output_stream.open(filename);
  this->Bulk_mesh_pt->output(output_stream, n_points);
  output_stream.close();
}

template<class ELEMENT>
void HeleShawChannelProblem<ELEMENT>::save_integral_to_file(
  ofstream& output_stream, string filename)
{
  output_stream.open(filename.c_str());
  double my_integral = 0.0;
  my_integral = Inlet_integral_data_pt->value(0);
  output_stream << "Inlet integral "
                << "= " << my_integral << endl;
  my_integral = Outlet_integral_data_pt->value(0);
  output_stream << "Outlet integral "
                << "= " << my_integral << endl;
  output_stream.close();
}

#endif
