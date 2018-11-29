#ifndef LIBMESH_FE
#define LIBMESH_FE

struct LibMesh_FE {
  // Finite Element type
  const FEType& fe_type;
  // Finite Element object.  Since the FEBase::build() member
  // dynamically creates memory we will store the object as a
  // std::unique_ptr<FEBase>.  This can be thought of as a pointer
  // that will clean up after itself.
  const int dim;
  // Inernal Finite Element object. Since the FEBase::build() member
  // dynamically creates memory we will store the object as a
  // std::unique_ptr<FEBase>. This can be thought of as a pointer
  // that will clean up after itself.
  std::unique_ptr<FEBase> fe;
  // Quadrature rule
  const QGauss& qrule;
  // Jacobian matrix times quadrature wheights
  const std::vector<Real>& JxW;
  // Basis functions evaluated on quadrature points
  const std::vector<std::vector<Real>> & phi;
  // Gradient of basis functions evaluated on quadrature points
  const std::vector<std::vector<RealGradient>> & dphi;
  // Quarature rule nodes
  const std::vector<Point> & qrule_points;
  // Normal vectors
  const std::vector<Point> & qface_normals;
  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  // Number of degrees of freedom. Set by reset() method.
  unsigned int n_dofs;

LibMesh_FE(int dim_, const FEType& fe_type_, const QGauss& qrule_) :
  dim(dim_),
  fe_type(fe_type_),
  fe(FEBase::build(dim, fe_type)),
  qrule(qrule_),
  JxW(fe->get_JxW()),
  phi(fe->get_phi()),
  dphi(fe->get_dphi()),
  qrule_points(fe->get_xyz()),
  qrule_normals(fe_elem_face->get_normals()),
  dof_indices(),
  n_dofs((unsigned int) -1);
  {
    fe.attach_quadrature_rule(&qrule);
  }

  // Reinit element from a DofMap
  void init_dofs(const Elem* element, const DofMap & dof_map,
		 bool reinit=true) {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices(element, dof_indices);
      n_dofs = dof_indices.size();

      if(reinit) {
	// Compute the element-specific data for the current
	// element.  This involves computing the location of the
	// quadrature points (q_point) and the shape functions
	// (phi, dphi) for the current element.
	fe->reinit(element);
      }
  }

};

#endif // LIBMESH_FE
