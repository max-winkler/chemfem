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
    class RaviartThomasElement : public VectorElement
    {
    public:
      RaviartThomasElement();

      Vector2D Value(int, double, double) const override;

      /// The identity for every basis function, they are all of the form x + const
      Matrix2D Gradient(int, double, double) const override;

      /// Midpoint of the edge the DOF belongs to
      chemfem::linalg::Coordinate NodalPoint(int) const override;

      /// The flux through the edge k, so an edge moment rather than a point value
      DofDescriptor Dof(int) const override;
    };
  };
};

#endif
