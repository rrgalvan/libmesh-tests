#ifndef DG_FACE_COUPLING
#define DG_FACE_COUPLING

class DG_FaceCoupling {
  // Coupling between an _element_ and its _neighbor_ through a common face.
  //
  // When we compute integrals in a face, e, of u*v, where
  //   u = restriction to e of _unknown_ (defined in element or neighbor),
  //   v = restriction to e of _test_ f. (defined in element or neighbor),
  // then results shall afect the FE global matrix in positions defined by
  // the local degrees of freedom of u and v in element and in neighbor.
  //
  // This class stores local integrals in element and neighbor of u and v,
  // making possible adding them to global matrix.
 public:
  typedef double (*TermFunction)(int,int);

  DG_FaceCoupling(const std::vector<libMesh::dof_id_type>& element_dof_indices,
		  const std::vector<libMesh::dof_id_type>& neighbor_dof_indices,
		  const std::vector<libMesh::Real>& JxW_face )
    : _element_dof_indices(element_dof_indices),
    _neighbor_dof_indices(neighbor_dof_indices),
    _n_element_dofs(element_dof_indices.size()),
    _n_neighbor_dofs(neighbor_dof_indices.size()),
    _JxW_face(JxW_face)
    /* _elem(elem), */
    /* _neighbor(elem->neighbor_ptr(side_id)), // Neighbor element */
    /* _side_id(side_id), */
      {}
  // Update local matrices from a function defining a variational formulation term
  // For example: this function could be similar to this way: f(i,j) -> u(i)*v(j) where
  //  u(i) = phi_face[i][qp], v(j) = phi_face[j][pq]
  template <typename ElementOrNeighborT1, typename ElementOrNeighborT2>
    void add_term(TermFunction) {
    throw std::logic_error("add_term<ElementOrNeighborT1, ElementOrNeighborT2>() must be specialized");
  }
  // Add local matrices to the global matrix of an implicit system
  void add_to_system(libMesh::LinearImplicitSystem & system);

 private:
  /* const Elem* _elem; */
  /* int _side_id; */
  const std::vector<libMesh::Real>& _JxW_face;
  const std::vector<libMesh::dof_id_type>& _element_dof_indices;
  const std::vector<libMesh::dof_id_type>& _neighbor_dof_indices;
  int _n_element_dofs;
  int _n_neighbor_dofs;

  // Data structures to contain the element and neighbor boundary matrix
  // contribution. This matrices will do ;the coupling between the dofs of
  // the element and those of his neighbors.
  // Ken: matrix coupling elem and neighbor dofs
  libMesh::DenseMatrix<libMesh::Number> _Kne;
  libMesh::DenseMatrix<libMesh::Number> _Ken;
  libMesh::DenseMatrix<libMesh::Number> _Kee;
  libMesh::DenseMatrix<libMesh::Number> _Knn;

  // Save term in a DenseMatrix
  inline void _add_term(TermFunction, libMesh::DenseMatrix<libMesh::Number>&, int, int);
};

// Add local matrices to the global matrix of an implicit system
void DG_FaceCoupling::add_to_system(libMesh::LinearImplicitSystem& system) {
  system.matrix->add_matrix(_Kee, _element_dof_indices);
  system.matrix->add_matrix(_Ken, _element_dof_indices, _neighbor_dof_indices);
  system.matrix->add_matrix(_Kne, _neighbor_dof_indices, _element_dof_indices);
  system.matrix->add_matrix(_Knn, _neighbor_dof_indices);
}

inline void DG_FaceCoupling::_add_term(TermFunction f,
				       libMesh::DenseMatrix<libMesh::Number>& M,
				       int n_dofs_test_function, int n_dofs_unknown) {
  assert(M.m() == n_dofs_test_function);
  assert(M.n() == n_dofs_unknown);
  for (int i_test=0; i_test<n_dofs_test_function; i_test++)
    for (int i_unknown=0; i_unknown<n_dofs_unknown; i_unknown++)
      M(i_test,i_unknown) += f(i_unknown, i_test);
}

class Element {};  // Utility class for template functions
class Neighbor {}; // Utility class for template functions

template <> inline void
DG_FaceCoupling::add_term<Element,Element>(typename DG_FaceCoupling::TermFunction f) {
  // Store the results of integration in a face of the
  // *element* test function i against the *element* test function j
  _add_term(f, _Kee, _n_element_dofs, _n_element_dofs);
}

template <> inline void
DG_FaceCoupling::add_term<Element,Neighbor>(typename DG_FaceCoupling::TermFunction f) {
  // Store the results of integration in a face of the
  // *element* test function i against the *neighbor* test function j
  _add_term(f, _Ken, _n_element_dofs, _n_neighbor_dofs);
}

template <> inline void
DG_FaceCoupling::add_term<Neighbor,Element>(typename DG_FaceCoupling::TermFunction f) {
  // Store the results of integration in a face of the
  // *neighbor* test function i against the *element* test function j
  _add_term(f, _Kne, _n_neighbor_dofs, _n_element_dofs);
}

template <> inline void
DG_FaceCoupling::add_term<Neighbor,Neighbor>(typename DG_FaceCoupling::TermFunction f) {
  // Store the results of integration in a face of the
  // *neighbor* test function i against the *neighbor* test function j
  _add_term(f, _Knn, _n_neighbor_dofs, _n_neighbor_dofs);
}

#endif // DG_FACE_COUPLING
