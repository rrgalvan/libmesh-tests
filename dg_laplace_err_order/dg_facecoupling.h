#ifndef DG_FACECOUPLING
#define DG_FACECOUPLING

#include "fe_wrapper.h"

typedef enum {PLUS, MINUS} CoupledElement;

// Coupling of two FE_Wrapper objets (denoted plus and minus) trough a
// common face
class DG_FaceCoupling
{
 public:
 DG_FaceCoupling(FE_Wrapper  const& fe_plus, FE_Wrapper const& fe_minus):
   _fe_plus(fe_plus), _fe_minus(fe_minus), _n_quad_points(_fe_plus.n_quad_points()) {
   assert(_fe_plus.n_quad_points() == _fe_minus.n_quad_points());
 }

  // Return one of the two coupled FE_Wrapper objects (PLUS or MINUS)
  template <CoupledElement>
  const FE_Wrapper& fe() const;

  // Add some value to one of the internal matrices
  template <CoupledElement Elem1, CoupledElement Elem2> void
    add_value_to_matrix(unsigned int i, unsigned int j, Number const& value);

  // Add a bilinear term to an specific Element x Elemenet matrix coupling
  template <CoupledElement Elem1, CoupledElement Elem2, class DG_BilinearTerm> void
    add_term(DG_BilinearTerm const& term, Number weight);

  // Add a bilinear term to internal matrix
  template <class DG_Term>
    void integrate_bilinear_term(DG_Term const& term);

  unsigned int n_quad_points() const { return _n_quad_points; }
  unsigned int qp() const { return _qp; }

 private:
  const FE_Wrapper& _fe_plus;
  const FE_Wrapper& _fe_minus;
  const unsigned int _n_quad_points;
  unsigned int _qp;
  DenseMatrix<Number> K_plus_plus;
  DenseMatrix<Number> K_plus_minus;
  DenseMatrix<Number> K_minus_plus;
  DenseMatrix<Number> K_minus_minus;
};

template <> inline FE_Wrapper const&
DG_FaceCoupling::fe<PLUS>() const { return _fe_plus; }

template <> inline FE_Wrapper const&
DG_FaceCoupling::fe<MINUS>() const { return _fe_minus; }

template <> void inline
DG_FaceCoupling::add_value_to_matrix<PLUS,PLUS>(unsigned int i, unsigned int j, Number const& value) {
  K_plus_plus(i,j) += value;
}
template <> void inline
DG_FaceCoupling::add_value_to_matrix<PLUS,MINUS>(unsigned int i, unsigned int j, Number const& value) {
  K_plus_minus(i,j) += value;
}
template <> void inline
DG_FaceCoupling::add_value_to_matrix<MINUS,PLUS>(unsigned int i, unsigned int j, Number const& value) {
  K_minus_plus(i,j) += value;
}
template <> void inline
DG_FaceCoupling::add_value_to_matrix<MINUS,MINUS>(unsigned int i, unsigned int j, Number const& value) {
  K_minus_minus(i,j) += value;
}

template <CoupledElement Elem1, CoupledElement Elem2, class DG_BilinearTerm> void
DG_FaceCoupling::add_term(DG_BilinearTerm const& term, Number weight)
{
  const int n_dofs_1 = fe<Elem1>().n_dofs();
  for (unsigned int i=0; i<n_dofs_1; i++)
    {
      const int n_dofs_2 = fe<Elem2>().n_dofs();
      for (unsigned int j=0; j<n_dofs_2; j++)
	{
	  add_value_to_matrix<Elem1,Elem2>(i, j,
					   weight * term.value<Elem1,Elem2>(i,j));
	}
    }
}

template <class DG_Term> void
DG_FaceCoupling::integrate_bilinear_term(DG_Term const& term)
{
  for (_qp=0; _qp<n_quad_points(); _qp++) // Integrate on face
    {
      const auto weight = fe<PLUS>().JxW[_qp];
      add_term<PLUS, PLUS, DG_Term>(term, weight);
      add_term<PLUS, MINUS,DG_Term>(term, weight);
      add_term<MINUS,PLUS, DG_Term>(term, weight);
      add_term<MINUS,MINUS,DG_Term>(term, weight);
    }
}

#endif // DG_FACECOUPLING

// Local Variables:
// mode: c++
// End:
