template<class ELEMENT>
class HeleShawVolumeConstraintElement
  : public LineVolumeConstraintBoundingElement,
    public virtual FaceGeometry<ELEMENT>
{
public:
  /// \short Contructor: Specify bulk element and index of face to which
  /// this face element is to be attached
  HeleShawVolumeConstraintElement(FiniteElement* const& element_pt,
                                  const int& face_index)
    : FaceGeometry<ELEMENT>(), LineVolumeConstraintBoundingElement()
  {
    // Attach the geometrical information to the element, by
    // making the face element from the bulk element
    element_pt->build_face_element(face_index, this);
    i_bubble = 0;
  }

  unsigned i_bubble;

  void fill_in_generic_residual_contribution_volume_constraint(
    Vector<double>& residuals)
  {
    // Add in the volume constraint term if required
    const int local_eqn = this->ptraded_local_eqn();
    if (local_eqn >= 0)
    {
      //        std::cout << "Adding vol contribution" << std::endl;

      double dV = 0;
      calculate_volume_contribution(dV);
      // Add to residual with sign chosen so that the volume is
      // positive when the elements bound the fluid
      residuals[local_eqn] += dV;

      //         std::cout << "Volume dof " << residuals[local_eqn] <<
      //         std::endl;
    }
  }

  void calculate_volume_contribution(double& dV)
  {
    // Find out how many nodes there are
    const unsigned n_node = this->nnode();

    // Set up memeory for the shape functions
    Shape psif(n_node);
    DShape dpsifds(n_node, 1);

    // Set the value of n_intpt
    const unsigned n_intpt = this->integral_pt()->nweight();

    // Storage for the local coordinate
    Vector<double> s(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the local coordinate at the integration point
      s[0] = this->integral_pt()->knot(ipt, 0);

      // Get the integral weight
      double W = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape function at the knot point
      this->dshape_local_at_knot(ipt, psif, dpsifds);

      // Get position and tangent vector
      Vector<double> interpolated_t1(2, 0.0);
      Vector<double> interpolated_x(2, 0.0);
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directional components
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_x[i] += this->nodal_position(l, i) * psif(l);
          interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
        }
      }

      // Calculate the length of the tangent Vector
      double tlength = interpolated_t1[0] * interpolated_t1[0] +
                       interpolated_t1[1] * interpolated_t1[1];

      // Set the Jacobian of the line element
      double J = sqrt(tlength);

      // Now calculate the normal Vector
      Vector<double> interpolated_n(2);
      this->outer_unit_normal(ipt, interpolated_n);

      // Assemble dot product
      double dot = 0.0;
      //                double y = interpolated_x[1];
      double b = 0;
      Problem_Parameter::channel_height_function(interpolated_x, b);
      for (unsigned k = 0; k < 1; k++)
      {
        dot += interpolated_x[k] * interpolated_n[k] * b;
      }

      dV += dot * W * J;
    }
  }


  /// Fill in contribution to residuals and Jacobian. This is specific
  /// to solid-based elements in which derivatives w.r.t. to nodal
  /// positions are evaluated by finite differencing
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    // Call the generic routine
    this->fill_in_generic_residual_contribution_volume_constraint(residuals);

    // Shape derivatives
    // Call the generic finite difference routine for the solid variables
    this->fill_in_jacobian_from_solid_position_by_fd(jacobian);
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    //      std::cout << "Volume" << std::endl;
    fill_in_contribution_to_jacobian(residuals, jacobian);
  }


  /// \short The "global" intrinsic coordinate of the element when
  /// viewed as part of a geometric object should be given by
  /// the FaceElement representation, by default
  double zeta_nodal(const unsigned& n,
                    const unsigned& k,
                    const unsigned& i) const
  {
    return FaceElement::zeta_nodal(n, k, i);
  }
};


template<class ELEMENT>
class HeleShawCoMXConstraintElement
  : public LineVolumeConstraintBoundingElement,
    public virtual FaceGeometry<ELEMENT>
{
public:
  /// \short Contructor: Specify bulk element and index of face to which
  /// this face element is to be attached
  HeleShawCoMXConstraintElement(FiniteElement* const& element_pt,
                                const int& face_index)
    : FaceGeometry<ELEMENT>(), LineVolumeConstraintBoundingElement()
  {
    // Attach the geometrical information to the element, by
    // making the face element from the bulk element
    element_pt->build_face_element(face_index, this);
  }

  void fill_in_generic_residual_contribution_volume_constraint(
    Vector<double>& residuals)
  {
    // Add in the volume constraint term if required
    const int local_eqn = this->ptraded_local_eqn();
    if (local_eqn >= 0)
    {
      //            std::cout << "CoM dof" << std::endl;
      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memeory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local coordinate
      Vector<double> s(1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the local coordinate at the integration point
        s[0] = this->integral_pt()->knot(ipt, 0);

        // Get the integral weight
        double W = this->integral_pt()->weight(ipt);

        // Call the derivatives of the shape function at the knot point
        this->dshape_local_at_knot(ipt, psif, dpsifds);

        // Get position and tangent vector
        Vector<double> interpolated_t1(2, 0.0);
        Vector<double> interpolated_x(2, 0.0);
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Calculate the length of the tangent Vector
        double tlength = interpolated_t1[0] * interpolated_t1[0] +
                         interpolated_t1[1] * interpolated_t1[1];

        // Set the Jacobian of the line element
        double J = sqrt(tlength);

        // Now calculate the normal Vector
        Vector<double> interpolated_n(2);
        this->outer_unit_normal(ipt, interpolated_n);
        double b = 0;
        Problem_Parameter::channel_height_function(interpolated_x, b);

        // Assemble dot product
        double dot = 0.0;
        for (unsigned k = 0; k < 1; k++)
        {
          dot += interpolated_x[k] * interpolated_x[k] * interpolated_n[k] * b;
        }

        // Add to residual with sign chosen so that the volume is
        // positive when the elements bound the fluid
        residuals[local_eqn] += dot * W * J * 0.5;
      }
    }
  }
  /// Fill in contribution to residuals and Jacobian. This is specific
  /// to solid-based elements in which derivatives w.r.t. to nodal
  /// positions are evaluated by finite differencing
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    // Call the generic routine
    this->fill_in_generic_residual_contribution_volume_constraint(residuals);

    // Shape derivatives
    // Call the generic finite difference routine for the solid variables
    this->fill_in_jacobian_from_solid_position_by_fd(jacobian);
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    //       std::cout << "CoM_X" << std::endl;
    fill_in_contribution_to_jacobian(residuals, jacobian);
  }


  /// \short The "global" intrinsic coordinate of the element when
  /// viewed as part of a geometric object should be given by
  /// the FaceElement representation, by default
  double zeta_nodal(const unsigned& n,
                    const unsigned& k,
                    const unsigned& i) const
  {
    return FaceElement::zeta_nodal(n, k, i);
  }
};

template<class ELEMENT>
class HeleShawCoMYConstraintElement
  : public LineVolumeConstraintBoundingElement,
    public virtual FaceGeometry<ELEMENT>
{
public:
  /// \short Contructor: Specify bulk element and index of face to which
  /// this face element is to be attached
  HeleShawCoMYConstraintElement(FiniteElement* const& element_pt,
                                const int& face_index)
    : FaceGeometry<ELEMENT>(), LineVolumeConstraintBoundingElement()
  {
    // Attach the geometrical information to the element, by
    // making the face element from the bulk element
    element_pt->build_face_element(face_index, this);
  }

  void fill_in_generic_residual_contribution_volume_constraint(
    Vector<double>& residuals)
  {
    // Add in the volume constraint term if required
    const int local_eqn = this->ptraded_local_eqn();
    if (local_eqn >= 0)
    {
      //            std::cout << "CoM dof" << std::endl;
      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memeory for the shape functions
      Shape psif(n_node);
      DShape dpsifds(n_node, 1);

      // Set the value of n_intpt
      const unsigned n_intpt = this->integral_pt()->nweight();

      // Storage for the local coordinate
      Vector<double> s(1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the local coordinate at the integration point
        s[0] = this->integral_pt()->knot(ipt, 0);

        // Get the integral weight
        double W = this->integral_pt()->weight(ipt);

        // Call the derivatives of the shape function at the knot point
        this->dshape_local_at_knot(ipt, psif, dpsifds);

        // Get position and tangent vector
        Vector<double> interpolated_t1(2, 0.0);
        Vector<double> interpolated_x(2, 0.0);
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directional components
          for (unsigned i = 0; i < 2; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psif(l);
            interpolated_t1[i] += this->nodal_position(l, i) * dpsifds(l, 0);
          }
        }

        // Calculate the length of the tangent Vector
        double tlength = interpolated_t1[0] * interpolated_t1[0] +
                         interpolated_t1[1] * interpolated_t1[1];

        // Set the Jacobian of the line element
        double J = sqrt(tlength);

        // Now calculate the normal Vector
        Vector<double> interpolated_n(2);
        this->outer_unit_normal(ipt, interpolated_n);

        double b = 0;
        Problem_Parameter::channel_height_function(interpolated_x, b);

        // Assemble dot product
        double dot =
          interpolated_x[0] * interpolated_x[1] * interpolated_n[0] * b;

        // Add to residual with sign chosen so that the volume is
        // positive when the elements bound the fluid
        residuals[local_eqn] += dot * W * J;
      }
    }
  }
  /// Fill in contribution to residuals and Jacobian. This is specific
  /// to solid-based elements in which derivatives w.r.t. to nodal
  /// positions are evaluated by finite differencing
  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
    // Call the generic routine
    this->fill_in_generic_residual_contribution_volume_constraint(residuals);

    // Shape derivatives
    // Call the generic finite difference routine for the solid variables
    this->fill_in_jacobian_from_solid_position_by_fd(jacobian);

    GeneralisedElement::fill_in_jacobian_from_external_by_fd(
      residuals, jacobian, false);
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    //       std::cout << "CoM_Y" << std::endl;
    fill_in_contribution_to_jacobian(residuals, jacobian);
  }


  /// \short The "global" intrinsic coordinate of the element when
  /// viewed as part of a geometric object should be given by
  /// the FaceElement representation, by default
  double zeta_nodal(const unsigned& n,
                    const unsigned& k,
                    const unsigned& i) const
  {
    return FaceElement::zeta_nodal(n, k, i);
  }
};

namespace oomph
{
  //==========================================================================
  /// We want to calculate three integral measures (could do area too):
  /// the bubble volume
  /// and the two components of the bubble centre of mass.
  /// The pressure is set by the volume constraint,
  /// and the frame speed is set by the x-component of the centre of mass.
  /// But the y-component of the centre of mass is a free variable.
  /// And this third equation is self-referential, requiring external data
  /// dependence in the Jacobian.
  //=========================================================================

  class SelfReferentialVolumeConstraintElement : public VolumeConstraintElement
  {
  public:
    /// \short Constructor: Pass pointer to target volume, pointer to Data
    /// item whose value specified by index_of_traded_pressure represents
    /// the "Pressure" value that "traded" for the volume contraint.
    /// The Data is stored as external Data for this element.
    SelfReferentialVolumeConstraintElement(
      double* prescribed_volume_pt,
      Data* p_traded_data_pt,
      const unsigned& index_of_traded_pressure)
      : VolumeConstraintElement(
          prescribed_volume_pt, p_traded_data_pt, index_of_traded_pressure)
    {
    }


    /// \short Fill in the residuals for the volume constraint
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      VolumeConstraintElement::fill_in_contribution_to_residuals(residuals);
    }

    /// \short Fill in the residuals and jacobian for the volume constraint
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      VolumeConstraintElement::fill_in_contribution_to_jacobian(residuals,
                                                                jacobian);
      GeneralisedElement::fill_in_jacobian_from_external_by_fd(
        residuals, jacobian, false);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      //       std::cout << "SelfReferential" << std::endl;
      VolumeConstraintElement::fill_in_contribution_to_jacobian(residuals,
                                                                jacobian);
      GeneralisedElement::fill_in_jacobian_from_external_by_fd(
        residuals, jacobian, false);
    }
  };


} // namespace oomph
