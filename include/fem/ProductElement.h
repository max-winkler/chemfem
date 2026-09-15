#ifndef _PRODUCT_ELEMENT_H_
#define _PRODUCT_ELEMENT_H_

#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    /**
     * The Cartesian product of a scalar element with itself, the usual way to a vector
     * valued space:
     *
     * \code
     *   LagrangeElement P2(2);
     *   ProductElement Velocity(P2, 2);     // [P2]^2
     * \endcode
     *
     * The basis functions are e_c phi_k, so each of them has one component that does not
     * vanish. Value, Gradient and Hessian give that component, Component(k) says which one
     * it is. Nothing is computed here, the work stays in the scalar element.
     *
     * Elements whose basis functions are genuinely vector valued, like Raviart-Thomas or
     * Nedelec, do not fit this pattern. They need their own class with their own mapping.
     */
    class ProductElement : public ScalarElement
    {
    public:
      ProductElement(const ScalarElement&, int components);

      double Value(int, double, double) const;

      Vector2D Gradient(int, double, double) const;

      Matrix2D Hessian(int, double, double) const;

      chemfem::linalg::Coordinate NodalPoint(int) const;

      /// The scalar element the product is built from
      const ScalarElement& ScalarPart() const;

    private:
      const ScalarElement& Scalar;
    };

  };
};

#endif
