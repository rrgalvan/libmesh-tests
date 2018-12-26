#ifndef DG_FACETERM
#define DG_FACETERM

#include "dg_facecoupling.h"

class DG_FaceTerm {
public:
  DG_FaceTerm(const DG_FaceCoupling& face_coupling) :
    _face(face_coupling)
  {}
  // void set_dg_face_coupling(DG_FaceCoupling *face_coupling) {
  //   _face=face_coupling;
  // }
protected:
  const DG_FaceCoupling& _face;
};

// template<CoupledElement Elem1, CoupledElement Elem2, typename ConcreteFaceTerm> void
// DG_FaceTerm<ConcreteFaceTerm>::integrate_on_coupling(DG_FaceCoupling& face_coupling) {
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

class SIP_BilinearForm: public DG_FaceTerm
{
public:
  SIP_BilinearForm(const DG_FaceCoupling& face_coupling,
		   double penalty, double h_elem):
    DG_FaceTerm(face_coupling),
    _penalty(penalty),
    _h_elem(h_elem) {}

  template <CoupledElement Elem1, CoupledElement Elem2>
  Number value(unsigned int i, unsigned int j);
private:
  const double _penalty;
  const double _h_elem;
};

template <> Number
SIP_BilinearForm::value<PLUS,PLUS>(unsigned int unknown_id, unsigned int test_id)
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<PLUS>().phi[unknown_id] [_face.qp()];
  auto const& du = _face.fe<PLUS>().dphi[unknown_id][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<PLUS>().phi[test_id] [_face.qp()];
  auto const& dv = _face.fe<PLUS>().dphi[test_id][_face.qp()]; // Gradient

  return -0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
}

template <> Number
SIP_BilinearForm::value<PLUS,MINUS>(unsigned int unknown_id, unsigned int test_id)
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<PLUS>().phi[unknown_id] [_face.qp()];
  auto const& du = _face.fe<PLUS>().dphi[unknown_id][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<MINUS>().phi[test_id] [_face.qp()];
  auto const& dv = _face.fe<MINUS>().dphi[test_id][_face.qp()]; // Gradient

  return 0.5 * (v*(du*normal) - u*(dv*normal)) - (_penalty/_h_elem)*u*v;
}

template <> Number
SIP_BilinearForm::value<MINUS,PLUS>(unsigned int unknown_id, unsigned int test_id)
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<MINUS>().phi[unknown_id] [_face.qp()];
  auto const& du = _face.fe<MINUS>().dphi[unknown_id][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<PLUS>().phi[test_id] [_face.qp()];
  auto const& dv = _face.fe<PLUS>().dphi[test_id][_face.qp()]; // Gradient

  return 0.5 * (-v*(du*normal) + u*(dv*normal)) - (_penalty/_h_elem)*u*v;
}

template <> Number
SIP_BilinearForm::value<MINUS,MINUS>(unsigned int unknown_id, unsigned int test_id)
{
  // Normal vector, considered in PLUS orientation
  auto const& normal = _face.fe<PLUS>().qrule_normals[_face.qp()];
  // Unknown
  auto const& u  = _face.fe<MINUS>().phi[unknown_id] [_face.qp()];
  auto const& du = _face.fe<MINUS>().dphi[unknown_id][_face.qp()]; // Gradient
  // Test function
  auto const& v  = _face.fe<MINUS>().phi[test_id] [_face.qp()];
  auto const& dv = _face.fe<MINUS>().dphi[test_id][_face.qp()]; // Gradient

  return 0.5 * (v*(du*normal) + u*(dv*normal)) + (_penalty/_h_elem)*u*v;
}

#endif // DG_FACETERM

// Local Variables:
// mode: c++
// End:
