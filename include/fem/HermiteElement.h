#ifndef _HERMITE_ELEMENT_H_
#define _HERMITE_ELEMENT_H_

#include "fem/Element.h"

namespace chemfem {
namespace fem {

/**
 * Cubic Hermite element. The local DOFs follow the blocks of DofManager: per vertex the
 * function value and the derivatives along the two edges meeting there, so 0,1,2 belong to
 * vertex 0, 3,4,5 to vertex 1 and 6,7,8 to vertex 2, and 9 is the value in the barycenter.
 * The derivatives are taken along the edge vectors, not along unit directions, which is what
 * makes the element affine equivalent and what the DOF transformation expects.
 */
class HermiteElement : public ScalarElement {
 public:
  HermiteElement();
  /**
   * Returns the function value of the ansatz function on the reference element.
   */
  double Value(int, double, double) const override;

  Vector2D Gradient(int, double, double) const override;

  Matrix2D Hessian(int, double, double) const override;

  chemfem::linalg::Coordinate NodalPoint(int) const override;

  DofDescriptor Dof(int) const override;
};

}  // namespace fem
}  // namespace chemfem

#endif
