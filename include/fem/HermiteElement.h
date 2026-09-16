#ifndef _HERMITE_ELEMENT_H_
#define _HERMITE_ELEMENT_H_

#include "fem/Element.h"

namespace chemfem {
namespace fem {

class HermiteElement : public ScalarElement {
 public:
  HermiteElement();
  /**
   * Returns the function value of the ansatz function on the reference element.
   */
  double Value(int, double, double) const;

  Vector2D Gradient(int, double, double) const;

  Matrix2D Hessian(int, double, double) const;

  chemfem::linalg::Coordinate NodalPoint(int) const;
};

}  // namespace fem
}  // namespace chemfem

#endif
