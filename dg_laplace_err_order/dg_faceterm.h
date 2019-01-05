#ifndef DG_FACETERM
#define DG_FACETERM

#include "dg_facecoupling.h"

class DG_Term {
public:
  DG_Term(const DG_FaceCoupling& face_coupling) :
    _face(face_coupling)
  {}
  // void set_dg_face_coupling(DG_FaceCoupling *face_coupling) {
  //   _face=face_coupling;
  // }
protected:
  const DG_FaceCoupling& _face;
};

// DG_Term<ConcreteTerm>::integrate_on_coupling(DG_FaceCoupling& face_coupling) {
//   const int n_dofs_1 = face.fe<Elem1>().n_dofs;
//   for (unsigned int i=0; i<n_dofs_1; i++)
//     {
//       const int n_dofs_2 = fe<Elem2>().n_dofs;
//       for (unsigned int j=0; j<n_dofs_2; j++)
// 	{
// 	  face_coupling.add_value_to_matrix<Elem1,Elem2>(i, j, value<Elem1,Elem2>(i,j));
//
//     }
// }

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
  Number inline value_plus_plus(unsigned int id_u, unsigned int id_test) const;
  Number inline value_plus_minus(unsigned int id_u, unsigned int id_test) const;
  Number inline value_minus_plus(unsigned int id_u, unsigned int id_test) const;
  Number inline value_minus_minus(unsigned int id_u, unsigned int id_test) const;

  // template <CoupledElement Elem1, CoupledElement Elem2>
  // Number value(unsigned int id_u, unsigned int id_test) const;
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


template <> Number inline
SIP_BilinearForm::value_plus_plus(unsigned int id_u, unsigned int id_test) const
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
SIP_BilinearForm::value_plus_minus(unsigned int id_u, unsigned int id_test) const
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

template <> Number
SIP_BilinearForm::value_minus_plus(unsigned int id_u, unsigned int id_test) const
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

template <> Number
SIP_BilinearForm::value_minus_minus(unsigned int id_u, unsigned int id_test) const
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

#endif // DG_FACETERM

// Local Variables:
// mode: c++
// End:
