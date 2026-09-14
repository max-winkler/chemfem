#ifndef _RAVIART_THOMAS_ELEMENT_H_
#define _RAVIART_THOMAS_ELEMENT_H_

#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    /**
     * Raviart-Thomas element of lowest order. Its DOF on an edge is the flux of the function
     * through that edge, so the normal component is continuous between two cells while the
     * tangential one jumps. The reference basis functions are
     *
     *   s_0 = (xi, eta-1),   s_1 = (xi, eta),   s_2 = (xi-1, eta),
     *
     * numbered like the edges, each with flux 1 through its own edge, flux 0 through the
     * other two and divergence 2. A cell maps them with the contravariant Piola transform,
     * so the scalar Value of the base class does not describe them.
     */
    class RaviartThomasElement : public Element
    {
    public:
      RaviartThomasElement();

      VectorRefValues VectorReference(int, double, double) const;

      /// Not available, the basis functions are vector valued
      double Value(int, double, double) const;

      /// Not available, the basis functions are vector valued
      Vector2D Gradient(int, double, double) const;

      /// Not available, the basis functions are vector valued
      Matrix2D Hessian(int, double, double) const;

      /// Midpoint of the edge the DOF belongs to
      chemfem::linalg::Coordinate NodalPoint(int) const;
    };
  };
};

#endif
