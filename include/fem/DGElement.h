#ifndef _DG_ELEMENT_H_
#define _DG_ELEMENT_H_

#include "fem/Element.h"

namespace chemfem{
  namespace fem{

    /**
     * Discontinuous Galerkin element: the polynomials of the given degree on each cell,
     * with no continuity requirement between the cells. All of its DOFs are therefore
     * interior ones. Only the degree 0 is implemented, the piecewise constants. They are
     * the pressure space of the lowest order Raviart-Thomas discretization, for which
     * div(RT_0) = P_0 holds exactly.
     */
    class DGElement : public ScalarElement
    {
    public:
      DGElement(int degree);

      double Value(int, double, double) const override;

      Vector2D Gradient(int, double, double) const override;

      Matrix2D Hessian(int, double, double) const override;

      chemfem::linalg::Coordinate NodalPoint(int) const override;
    };
  };
};

#endif
