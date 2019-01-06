#ifndef DG_FACECOUPLING
#define DG_FACECOUPLING

#include "fe_wrapper.h"

typedef enum {PLUS, MINUS} CoupledElement;

// Coupling of two FE_Wrapper objets (denoted plus and minus) trough a
// common face
class DG_FaceCoupling
{
 public:
 DG_FaceCoupling(FE_Wrapper const& fe_plus, FE_Wrapper const& fe_minus):
   _fe_plus(fe_plus), _fe_minus(fe_minus), _n_quad_points(_fe_plus.n_quad_points())
  {
    assert(_fe_plus.n_quad_points() == _fe_minus.n_quad_points());
  }

  // Return one of the two coupled FE_Wrapper objects (PLUS or MINUS)
  template <CoupledElement>
  const FE_Wrapper& fe() const;

  // // Add a bilinear term to an specific Element x Elemenet matrix coupling
  // template <CoupledElement Elem1, CoupledElement Elem2, class DG_BilinearTerm> void
  //   add_term(DG_BilinearTerm const& term, Number weight);

  // // Add a bilinear term to internal matrix
  // template <class DG_Term>
  //   void integrate_bilinear_term(DG_Term const& term);

  unsigned int n_quad_points() const { return _n_quad_points; }
  unsigned int qp() const { return _qp; }

  // Save local matrices and RHS to global values in a System
  void add_to_system(LinearImplicitSystem& sys);

 private:
  const FE_Wrapper& _fe_plus;
  const FE_Wrapper& _fe_minus;
  const unsigned int _n_quad_points;
  unsigned int _qp;
};

template <> inline FE_Wrapper const&
DG_FaceCoupling::fe<PLUS>() const { return _fe_plus; }

template <> inline FE_Wrapper const&
DG_FaceCoupling::fe<MINUS>() const { return _fe_minus; }

//----------------------------------------------------------------------

class DG_Term
{
public:
  DG_Term(const DG_FaceCoupling& face_coupling) :
    _face(face_coupling) {}
  DG_FaceCoupling const& face() const { return _face; };
protected:
  const DG_FaceCoupling& _face;
};

class SIP_BilinearForm: public DG_Term
{
public:
  SIP_BilinearForm(const DG_FaceCoupling& face_coupling,
		   double penalty, double h_elem):
    DG_Term(face_coupling),
    _penalty(penalty),
    _h_elem(h_elem) {}

  // // Integrate term on a face
  // void integrate(DG_FaceCoupling*);

  // Evalate bilinear form (at current qudrature point) for the
  // unknown basis function id_u and the test basis function id_test
  Number inline value_plus_plus(unsigned int qp, unsigned int id_u, unsigned int id_test) const;
  Number inline value_plus_minus(unsigned int qp, unsigned int id_u, unsigned int id_test) const;
  Number inline value_minus_plus(unsigned int qp, unsigned int id_u, unsigned int id_test) const;
  Number inline value_minus_minus(unsigned int qp, unsigned int id_u, unsigned int id_test) const;

  template <CoupledElement Elem1, CoupledElement Elem2>
  Number value(unsigned int id_u, unsigned int id_test) const;
private:
  const double _penalty;
  const double _h_elem;
};

template <> Number inline
SIP_BilinearForm::value<PLUS,PLUS>(unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<PLUS>().phi[id_u] [_face.qp()];
  auto const& du = _face.fe<PLUS>().dphi[id_u][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<PLUS>().phi[id_test] [_face.qp()];
  auto const& dv = _face.fe<PLUS>().dphi[id_test][_face.qp()]; // Gradient

  return -0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
}

template <> Number inline
SIP_BilinearForm::value<PLUS,MINUS>(unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<PLUS>().phi[id_u] [_face.qp()];
  auto const& du = _face.fe<PLUS>().dphi[id_u][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<MINUS>().phi[id_test] [_face.qp()];
  auto const& dv = _face.fe<MINUS>().dphi[id_test][_face.qp()]; // Gradient

  return 0.5 * (v*(du*normal) - u*(dv*normal)) - (_penalty/_h_elem)*u*v;
}

template <> Number inline
SIP_BilinearForm::value<MINUS,PLUS>(unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<MINUS>().phi[id_u] [_face.qp()];
  auto const& du = _face.fe<MINUS>().dphi[id_u][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<PLUS>().phi[id_test] [_face.qp()];
  auto const& dv = _face.fe<PLUS>().dphi[id_test][_face.qp()]; // Gradient

  return 0.5 * (-v*(du*normal) + u*(dv*normal)) - (_penalty/_h_elem)*u*v;
}

template <> Number inline
SIP_BilinearForm::value<MINUS,MINUS>(unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<MINUS>().phi[id_u] [_face.qp()];
  auto const& du = _face.fe<MINUS>().dphi[id_u][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<MINUS>().phi[id_test] [_face.qp()];
  auto const& dv = _face.fe<MINUS>().dphi[id_test][_face.qp()]; // Gradient

  return 0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
}

//----------------------------------------------------------------------


Number inline
SIP_BilinearForm::value_plus_plus(unsigned int qp, unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[qp];
  // Unknown
  auto const& u  = _face.fe<PLUS>().phi[id_u] [qp];
  auto const& du = _face.fe<PLUS>().dphi[id_u][qp]; // Gradient
  // Test function
  auto const& v  = _face.fe<PLUS>().phi[id_test] [qp];
  auto const& dv = _face.fe<PLUS>().dphi[id_test][qp]; // Gradient

  return -0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
}

Number inline
SIP_BilinearForm::value_plus_minus(unsigned int qp, unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[qp];
  // Unknown
  auto const& u  = _face.fe<PLUS>().phi[id_u] [qp];
  auto const& du = _face.fe<PLUS>().dphi[id_u][qp]; // Gradient
  // Test function
  auto const& v  = _face.fe<MINUS>().phi[id_test] [qp];
  auto const& dv = _face.fe<MINUS>().dphi[id_test][qp]; // Gradient

  return 0.5 * (v*(du*normal) - u*(dv*normal)) - (_penalty/_h_elem)*u*v;
}

Number inline
SIP_BilinearForm::value_minus_plus(unsigned int qp, unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[qp];
  // Unknown
  auto const& u  = _face.fe<MINUS>().phi[id_u] [qp];
  auto const& du = _face.fe<MINUS>().dphi[id_u][qp]; // Gradient
  // Test function
  auto const& v  = _face.fe<PLUS>().phi[id_test] [qp];
  auto const& dv = _face.fe<PLUS>().dphi[id_test][qp]; // Gradient

  return 0.5 * (-v*(du*normal) + u*(dv*normal)) - (_penalty/_h_elem)*u*v;
}

Number inline
SIP_BilinearForm::value_minus_minus(unsigned int qp, unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[qp];
  // Unknown
  auto const& u  = _face.fe<MINUS>().phi[id_u] [qp];
  auto const& du = _face.fe<MINUS>().dphi[id_u][qp]; // Gradient
  // Test function
  auto const& v  = _face.fe<MINUS>().phi[id_test] [qp];
  auto const& dv = _face.fe<MINUS>().dphi[id_test][qp]; // Gradient

  return 0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
}

// //----------------------------------------------------------------------

// class DG_Term
// {
// public:
//   DG_Term(const DG_FaceCoupling& face_coupling) :
//     _face(face_coupling) {}
//   DG_FaceCoupling const& face() const { return _face; };
// protected:
//   const DG_FaceCoupling& _face;
// };

// class SIP_BilinearForm: public DG_Term
// {
// public:
//   SIP_BilinearForm(const DG_FaceCoupling& face_coupling,
// 		   double penalty, double h_elem):
//     DG_Term(face_coupling),
//     _penalty(penalty),
//     _h_elem(h_elem) {}


//   template <CoupledElement Elem1, CoupledElement Elem2>
//   Number value(unsigned int i, unsigned int j) const;
// private:
//   const double _penalty;
//   const double _h_elem;
// };

// // void SIP_BilinearForm::integrate(DG_FaceCoupling* face)
// // {
// //   for (unsigned int qp=0; qp<_face.n_quad_points(); qp++) // Integrate on face
// //     {
// //       const auto weight = _face.fe<PLUS>().JxW[qp];
// //       face->add_term<PLUS, PLUS, SIP_BilinearForm>(*this, weight);
// //       face->add_term<PLUS, MINUS,SIP_BilinearForm>(*this, weight);
// //       face->add_term<MINUS,PLUS, SIP_BilinearForm>(*this, weight);
// //       face->add_term<MINUS,MINUS,SIP_BilinearForm>(*this, weight);
// //     }
// // }

// template <> Number
// SIP_BilinearForm::value<PLUS,PLUS>(unsigned int unknown_id, unsigned int test_id) const
// {
//   // Normal vector, considered in PLUS orientation
//   auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
//   // Unknown
//   auto const& u  = _face.fe<PLUS>().phi[unknown_id] [_face.qp()];
//   auto const& du = _face.fe<PLUS>().dphi[unknown_id][_face.qp()]; // Gradient
//   // Test function
//   auto const& v  = _face.fe<PLUS>().phi[test_id] [_face.qp()];
//   auto const& dv = _face.fe<PLUS>().dphi[test_id][_face.qp()]; // Gradient

//   return -0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
// }

// template <> Number
// SIP_BilinearForm::value<PLUS,MINUS>(unsigned int unknown_id, unsigned int test_id) const
// {
//   // Normal vector, considered in PLUS orientation
//   auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
//   // Unknown
//   auto const& u  = _face.fe<PLUS>().phi[unknown_id] [_face.qp()];
//   auto const& du = _face.fe<PLUS>().dphi[unknown_id][_face.qp()]; // Gradient
//   // Test function
//   auto const& v  = _face.fe<MINUS>().phi[test_id] [_face.qp()];
//   auto const& dv = _face.fe<MINUS>().dphi[test_id][_face.qp()]; // Gradient

//   return 0.5 * (v*(du*normal) - u*(dv*normal)) - (_penalty/_h_elem)*u*v;
// }

// template <> Number
// SIP_BilinearForm::value<MINUS,PLUS>(unsigned int unknown_id, unsigned int test_id) const
// {
//   // Normal vector, considered in PLUS orientation
//   auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
//   // Unknown
//   auto const& u  = _face.fe<MINUS>().phi[unknown_id] [_face.qp()];
//   auto const& du = _face.fe<MINUS>().dphi[unknown_id][_face.qp()]; // Gradient
//   // Test function
//   auto const& v  = _face.fe<PLUS>().phi[test_id] [_face.qp()];
//   auto const& dv = _face.fe<PLUS>().dphi[test_id][_face.qp()]; // Gradient

//   return 0.5 * (-v*(du*normal) + u*(dv*normal)) - (_penalty/_h_elem)*u*v;
// }

// template <> Number
// SIP_BilinearForm::value<MINUS,MINUS>(unsigned int unknown_id, unsigned int test_id) const
// {
//   // Normal vector, considered in PLUS orientation
//   auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
//   // Unknown
//   auto const& u  = _face.fe<MINUS>().phi[unknown_id] [_face.qp()];
//   auto const& du = _face.fe<MINUS>().dphi[unknown_id][_face.qp()]; // Gradient
//   // Test function
//   auto const& v  = _face.fe<MINUS>().phi[test_id] [_face.qp()];
//   auto const& dv = _face.fe<MINUS>().dphi[test_id][_face.qp()]; // Gradient

//   return 0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
// }

//----------------------------------------------------------------------

template <CoupledElement Elem1, CoupledElement Elem2>
class face_select {};

template <class Term>
class FaceIntegrator
{
public:
  FaceIntegrator(Term const& t): _face(t.face()), _term(t)
  {
    _K_plus_plus.resize   ( _face.fe<PLUS>().n_dofs(), _face.fe<PLUS>().n_dofs() );
    _K_plus_minus.resize  ( _face.fe<PLUS>().n_dofs(), _face.fe<MINUS>().n_dofs() );
    _K_minus_plus.resize  ( _face.fe<MINUS>().n_dofs(), _face.fe<PLUS>().n_dofs() );
    _K_minus_minus.resize ( _face.fe<MINUS>().n_dofs(), _face.fe<MINUS>().n_dofs() );
  }
  FaceIntegrator(DG_FaceCoupling const& f, Term const& t): _face(f), _term(t)
  {
    _K_plus_plus.resize   ( _face.fe<PLUS>().n_dofs(), _face.fe<PLUS>().n_dofs() );
    _K_plus_minus.resize  ( _face.fe<PLUS>().n_dofs(), _face.fe<MINUS>().n_dofs() );
    _K_minus_plus.resize  ( _face.fe<MINUS>().n_dofs(), _face.fe<PLUS>().n_dofs() );
    _K_minus_minus.resize ( _face.fe<MINUS>().n_dofs(), _face.fe<MINUS>().n_dofs() );
  }

  DenseMatrix<Number>& get_K_plus_plus()   { return _K_plus_plus; }
  DenseMatrix<Number>& get_K_plus_minus()  { return _K_plus_minus; }
  DenseMatrix<Number>& get_K_minus_plus()  { return _K_minus_plus; }
  DenseMatrix<Number>& get_K_minus_minus() { return _K_minus_minus; }

  void integrate();
  void save_to_system(LinearImplicitSystem& system);

private:
  DenseMatrix<Number> _K_plus_plus;
  DenseMatrix<Number> _K_plus_minus;
  DenseMatrix<Number> _K_minus_plus;
  DenseMatrix<Number> _K_minus_minus;
  DG_FaceCoupling const& _face;
  Term const& _term;
  unsigned int _qp;

  template <CoupledElement Elem1, CoupledElement Elem2> void
    add_values_to_matrix(Number const& weight);

};

template <class Term>
template <CoupledElement Elem1, CoupledElement Elem2> void
FaceIntegrator<Term>::add_values_to_matrix(Number const& weight)
{
  const unsigned int n_dofs_1 = _face.fe<Elem1>().n_dofs();
  const unsigned int n_dofs_2 = _face.fe<Elem2>().n_dofs();
  for (unsigned int i=0; i<n_dofs_1; i++)
    {
      for (unsigned int j=0; j<n_dofs_2; j++)
	{
	  // Number value = _term.value<Elem1,Elem2>(i,j);
	  Number value = face_select<Elem1,Elem2>::evaluate(_term,_qp,i,j);
	  DenseMatrix<Number>& K = face_select<Elem1,Elem2>::matrix(this);
	  K(i,j) += weight*value;
	}
    }
}

template <class BilinearForm> void
FaceIntegrator<BilinearForm>::integrate()
{
  std::cout << "### 3000" << " n_qp=" << _face.n_quad_points() << std::endl;
  for (_qp=0; _qp<_face.n_quad_points(); _qp++) // Integrate on face
    {
      std::cout << "### 3001" << " qp=" << _qp << std::endl;
      const auto weight = _face.fe<PLUS>().JxW[_qp];
      std::cout << "### 3002" << std::endl;
      add_values_to_matrix<PLUS, PLUS> (weight);
      std::cout << "### 3003" << std::endl;
      add_values_to_matrix<PLUS, MINUS>(weight);
      add_values_to_matrix<MINUS,PLUS> (weight);
      add_values_to_matrix<MINUS,MINUS>(weight);
    }
}

template <class BilinearForm> void
FaceIntegrator<BilinearForm>::save_to_system(LinearImplicitSystem& system)
{
  const auto& dofs_plus  = _face.fe<PLUS>().dof_indices;
  const auto& dofs_minus = _face.fe<MINUS>().dof_indices;
  system.matrix->add_matrix(_K_plus_plus, dofs_plus);
  system.matrix->add_matrix(_K_plus_minus, dofs_plus, dofs_minus);
  system.matrix->add_matrix(_K_minus_plus, dofs_minus, dofs_plus);
  system.matrix->add_matrix(_K_minus_minus, dofs_minus);
}

// ----------------------------------------------------

template <>
struct face_select<PLUS,PLUS> {
  template <class Term>
  static DenseMatrix<Number>& matrix(FaceIntegrator<Term>* fi)
  {
    return fi->get_K_plus_plus();
  }
  template <class Term>
  static Number evaluate(const Term& term, unsigned int qp, unsigned int i, unsigned int j)
  {
    return term.value_plus_plus(qp,i,j);
  }

};

template <>
struct face_select<PLUS,MINUS> {
  template <class Term>
  static DenseMatrix<Number>& matrix(FaceIntegrator<Term>* fi)
  {
    return fi->get_K_plus_minus();
  }
  template <class Term>
  static Number evaluate(const Term& term, unsigned int qp, unsigned int i, unsigned int j)
  {
    return term.value_plus_minus(qp, i,j);
  }
};

template <>
struct face_select<MINUS,PLUS> {
  template <class Term>
  static DenseMatrix<Number>& matrix(FaceIntegrator<Term>* fi)
  {
    return fi->get_K_minus_plus();
  }
  template <class Term>
  static Number evaluate(const Term& term, unsigned int qp, unsigned int i, unsigned int j)
  {
    return term.value_minus_plus(qp,i,j);
  }
};

template <>
struct face_select<MINUS,MINUS> {
  template <class Term>
  static DenseMatrix<Number>& matrix(FaceIntegrator<Term>* fi)
  {
    return fi->get_K_minus_minus();
  }
  template <class Term>
  static Number evaluate(const Term& term, unsigned int qp, unsigned int i, unsigned int j)
  {
    return term.value_minus_minus(qp,i,j);
  }
};

#endif // DG_FACECOUPLING


// Local Variables:
// mode: c++
// End:
