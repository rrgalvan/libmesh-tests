#ifndef FE_WRAPPER
#define FE_WRAPPER

#include "libmesh/libmesh.h"
using namespace libMesh;

//
// Class wrapping libmesh objects: finite element, quadrature rule, etc.
//
// The method init_dofs computes this->dof_indices the degrees of
// freedom indices or the current element, defining location of
// matrix and rhs contributions.
//
struct FE_Wrapper {
  // Inernal Finite Element object. Since the FEBase::build() member
  // dynamically creates memory we will store the object as a
  // std::unique_ptr<FEBase>. This can be thought of as a pointer
  // that will clean up after itself.
  std::unique_ptr<FEBase> fe;
  // Quadrature rule
  QGauss* qrule;
  // Jacobian matrix times quadrature wheights
  const std::vector<Real>& JxW;
  // Basis functions evaluated on quadrature points
  const std::vector<std::vector<Real>> & phi;
  // Gradient of basis functions evaluated on quadrature points
  const std::vector<std::vector<RealGradient>> & dphi;
  // Quarature rule nodes
  const std::vector<Point> & qrule_points;
  // Normal vectors
  const std::vector<Point> & qrule_normals;
  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  // Number of degrees of freedom. Set by reset() method.
  unsigned int _n_dofs;

FE_Wrapper(std::unique_ptr<FEBase>&& fe_, QGauss* qrule_):
  fe(std::move(fe_)),
  qrule(qrule_),
  JxW(fe->get_JxW()),
  phi(fe->get_phi()),
  dphi(fe->get_dphi()),
  qrule_points(fe->get_xyz()),
  qrule_normals(fe->get_normals()),
  dof_indices(),
  _n_dofs((unsigned int) -1)//,
  // n_quad_points(qrule_->n_points())
  {
    fe->attach_quadrature_rule(qrule_);
  }

  // Reinit element degrees of freedom from a DofMap
  void init_dofs(const Elem* element, const DofMap & dof_map)
  {
      // Compute (in this->dof_indices) the degrees of freedom indices
      // for the current element.  These define where in the global
      // matrix and right-hand-side this element will contribute to.
      dof_map.dof_indices(element, dof_indices);
      _n_dofs = dof_indices.size();
  }

  // Reinit element degrees of freedom from a given vector of dof indices
  void init_dofs(const std::vector<dof_id_type>& dof_indices)
  {
      this->dof_indices = dof_indices;
      _n_dofs = dof_indices.size();
  }

  unsigned int n_dofs() const { return _n_dofs; };
  unsigned int n_quad_points() const { return qrule->n_points(); }

};

#endif // FE_WRAPPER

// Local Variables:
// mode: c++
// End:
