#pragma once

#include "generic.h"

namespace gp {
void ode(const double &t, const double &u, double &dudt) { dudt = -u; };
} // namespace gp

class ODEProblem : public Problem {
public:
  ODEProblem() {
    // Create and add time stepper
    TimeStepper *time_stepper_pt = new BDF<2>;
    this->add_time_stepper_pt(time_stepper_pt);

    // Create ODE element
    ode_element_pt = new MyODEElement(time_stepper_pt, &gp::ode);

    // Create mesh
    this->mesh_pt() = new Mesh;
    this->mesh_pt()->add_element_pt(ode_element_pt);

    // Setup equation numbering
    this->assign_eqn_numbers();
    cout << "Number of equations: " << this->ndof() << endl;
  }

  ~ODEProblem(){};

  void set_initial_condition() {
    unsigned element_index = 0;
    unsigned data_index = 0;
    unsigned value_index = 0;
    unsigned time_level = 0;
    double x0 = 2.0;

    this->mesh_pt()
        ->element_pt(element_index)
        ->internal_data_pt(data_index)
        ->set_value(time_level, value_index, x0);

    time_level = 1;
    x0 = 1.0;
    this->mesh_pt()
        ->element_pt(element_index)
        ->internal_data_pt(data_index)
        ->set_value(time_level, value_index, x0);
  }

  void iterate_timestepper(const double &t_step, const double &t_final,
                           DocInfo &doc_info) {
    double t_duration = t_final - this->time_stepper_pt()->time();
    unsigned n_timestep = ceil(t_duration / t_step);
    for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++) {
      this->unsteady_newton_solve(t_step);
      this->doc_solution(doc_info);
    }
  }

  void doc_solution(DocInfo &doc_info) {

    ofstream output_stream;
    string filename =
        doc_info.directory() + "/soln" + to_string(doc_info.number()) + ".dat";

    output_stream.open(filename);
    unsigned element_index = 0;
    dynamic_cast<MyODEElement *>(this->mesh_pt()->element_pt(element_index))
        ->output(output_stream);
    output_stream.close();

    doc_info.number()++;
  }

private:
  MyODEElement *ode_element_pt;
};
