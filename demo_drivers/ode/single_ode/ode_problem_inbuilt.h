#pragma once

#include "generic.h"
#include "ode.h"

namespace problem_functor
{
  class ODE : public SolutionFunctorBase
  {
  public:
    ODE() {}

    virtual ~ODE() {}

    /// Exact solution -- Not necessarily known, used for setting the initial
    /// values and initial past values
    Vector<double> operator()(const double& t, const Vector<double>& x) const
    {
      Vector<double> output(1);

      output[0] = std::exp(-t);

      return output;
    }

    /// First-order ODE right hand side
    Vector<double> derivative(const double& t,
                              const Vector<double>& x,
                              const Vector<double>& u) const
    {
      Vector<double> output(1);

      output[0] = -u[0];

      return output;
    }
  };
} // namespace problem_functor

class ODEProblem : public Problem
{
public:
  ODEProblem()
  {
    // Create and add time stepper
    TimeStepper* time_stepper_pt = new BDF<2>;
    this->add_time_stepper_pt(time_stepper_pt);

    // Create ODE element
    SolutionFunctorBase* exact_solution_pt = new problem_functor::ODE;
    ode_element_pt = new ODEElement(time_stepper_pt, exact_solution_pt);

    // Create mesh
    this->mesh_pt() = new Mesh;
    this->mesh_pt()->add_element_pt(ode_element_pt);

    // Setup equation numbering
    this->assign_eqn_numbers();
    cout << "Number of equations: " << this->ndof() << endl;
  }

  ~ODEProblem(){};

  void set_initial_condition()
  {
    unsigned element_index = 0;
    unsigned data_index = 0;
    unsigned value_index = 0;

    unsigned n_tstorage = this->time_stepper_pt()->ntstorage();
    for (unsigned time_level = 0; time_level < n_tstorage; time_level++)
    {
      double t = this->time_stepper_pt()->time(time_level);
      Vector<double> dummy_x(1);
      Vector<double> u = problem_functor::ODE(t, dummy_x);
      this->mesh_pt()
        ->element_pt(element_index)
        ->internal_data_pt(data_index)
        ->set_value(time_level, value_index, u);
    }
  }

  void iterate_timestepper(const double& t_step,
                           const double& t_final,
                           DocInfo& doc_info)
  {
    double t_duration = t_final - this->time_stepper_pt()->time();
    unsigned n_timestep = ceil(t_duration / t_step);
    for (unsigned i_timestep = 0; i_timestep < n_timestep; i_timestep++)
    {
      this->unsteady_newton_solve(t_step);
      this->doc_solution(doc_info);
    }
  }

  void doc_solution(DocInfo& doc_info)
  {
    ofstream output_stream;
    string filename =
      doc_info.directory() + "/soln" + to_string(doc_info.number()) + ".dat";

    output_stream.open(filename);
    unsigned element_index = 0;
    // dynamic_cast<MyODEElement *>(this->mesh_pt()->element_pt(element_index))
    //    ->output(output_stream);

    unsigned data_index = 0;
    double t = this->mesh_pt()
                 ->element_pt(element_index)
                 ->internal_data_pt(data_index)
                 ->time_stepper_pt()
                 ->time();
    unsigned value_index = 0;
    double x = this->mesh_pt()
                 ->element_pt(element_index)
                 ->internal_data_pt(data_index)
                 ->value(value_index);
    output_stream << "t: " << t << ", x: " << x << endl;
    output_stream.close();

    doc_info.number()++;
  }

private:
  // MyODEElement *ode_element_pt;
  ODEElement* ode_element_pt;
};
