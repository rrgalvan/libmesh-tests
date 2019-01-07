#ifndef DG_FACECOUPLING
#define DG_FACECOUPLING

#include "fe_wrapper.h"

typedef enum {PLUS, MINUS} CoupledElement;

template <CoupledElement> Number inline jump(Number x);
template <> Number inline jump<PLUS>(Number x) { return x; }
template <> Number inline jump<MINUS>(Number x) { return -x; }

Number inline mean(Number x) { return 0.5*x; }

//----------------------------------------------------------------------
// Coupling of two FE_Wrapper objets (denoted plus and minus) trough a
// common face
class DG_FaceCoupling
//----------------------------------------------------------------------
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
//----------------------------------------------------------------------
{
public:
  DG_Term(const DG_FaceCoupling& face_coupling) :
    _face(face_coupling) {}
  DG_FaceCoupling const& face() const { return _face; };
protected:
  const DG_FaceCoupling& _face;
};

//----------------------------------------------------------------------
class SIP_BilinearForm: public DG_Term
//----------------------------------------------------------------------
{
public:
  SIP_BilinearForm(const DG_FaceCoupling& face_coupling,
		   double penalty, double h_elem):
    DG_Term(face_coupling),
    _penalty(penalty),
    _h_elem(h_elem) {}

  // Evalate bilinear form (at current qudrature point) for the
  // unknown basis function id_u and the test basis function id_test
  Number inline value_plus_plus(unsigned int qp, unsigned int id_u, unsigned int id_test) const {
    return _value<PLUS,PLUS>(qp,id_u,id_test); }
  Number inline value_plus_minus(unsigned int qp, unsigned int id_u, unsigned int id_test) const {
    return _value<PLUS,MINUS>(qp,id_u,id_test); }
  Number inline value_minus_plus(unsigned int qp, unsigned int id_u, unsigned int id_test) const {
    return _value<MINUS,PLUS>(qp,id_u,id_test); }
  Number inline value_minus_minus(unsigned int qp, unsigned int id_u, unsigned int id_test) const {
    return _value<MINUS,MINUS>(qp,id_u,id_test); }

private:
  const double _penalty;
  const double _h_elem;

  template <CoupledElement Elem1, CoupledElement Elem2>
  Number _value(unsigned int qp, unsigned int id_u, unsigned int id_test) const;
};

template <CoupledElement Elem1, CoupledElement Elem2> inline
Number SIP_BilinearForm::_value(unsigned int qp, unsigned int id_u, unsigned int id_test) const
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[qp];
  // Unknown
  auto const& u  = _face.fe<Elem1>().phi[id_u] [qp];
  auto const& du = _face.fe<Elem1>().dphi[id_u][qp]; // Gradient
  // Test function
  auto const& v  = _face.fe<Elem2>().phi[id_test] [qp];
  auto const& dv = _face.fe<Elem2>().dphi[id_test][qp]; // Gradient

  return - ( mean(du*normal)*jump<Elem2>(v) + jump<Elem1>(u)*mean(dv*normal) )
    + (_penalty/_h_elem) * jump<Elem1>(u)*jump<Elem2>(v);
}


//----------------------------------------------------------------------
template <CoupledElement Elem1, CoupledElement Elem2>
class face_select {};
//----------------------------------------------------------------------

//----------------------------------------------------------------------
template <class Term>
class FaceIntegrator
//----------------------------------------------------------------------
{
public:
  FaceIntegrator(Term const& t): _face(t.face()), _term(t)
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
    _integrate_to_matrix(Number const& weight);

};

template <class Term>
template <CoupledElement Elem1, CoupledElement Elem2> void
FaceIntegrator<Term>::_integrate_to_matrix(Number const& weight)
{
  const unsigned int n_dofs_1 = _face.fe<Elem1>().n_dofs();
  const unsigned int n_dofs_2 = _face.fe<Elem2>().n_dofs();
  for (unsigned int i=0; i<n_dofs_1; i++)
    {
      for (unsigned int j=0; j<n_dofs_2; j++)
	{
	  // Number value = _term.value<Elem1,Elem2>(i,j);
	  Number value = face_select<Elem1,Elem2>::evaluate(_term, _qp, i, j);
	  DenseMatrix<Number>& K = face_select<Elem1,Elem2>::matrix(this);
	  K(i,j) += weight*value;
	}
    }
}

template <class BilinearForm> void
FaceIntegrator<BilinearForm>::integrate()
{
  for (_qp=0; _qp<_face.n_quad_points(); _qp++) // Integrate on face
    {
      const auto weight = _face.fe<PLUS>().JxW[_qp];
      _integrate_to_matrix<PLUS, PLUS> (weight);
      _integrate_to_matrix<PLUS, MINUS>(weight);
      _integrate_to_matrix<MINUS,PLUS> (weight);
      _integrate_to_matrix<MINUS,MINUS>(weight);
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
