#ifndef _CROUZEIX_RAVIART_ELEMENT_H_
#define _CROUZEIX_RAVIART_ELEMENT_H_

#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    /**
     * Nonconforming P1 element of Crouzeix and Raviart. Its DOFs are the function values
     * in the edge midpoints, the local DOF k belongs to the local edge k.
     */
    class CrouzeixRaviartElement : public ScalarElement
    {
    public:
      CrouzeixRaviartElement();

      double Value(int, double, double) const override;

      Vector2D Gradient(int, double, double) const override;

      Matrix2D Hessian(int, double, double) const override;

      chemfem::linalg::Coordinate NodalPoint(int) const override;
    };
  };
};

#endif
