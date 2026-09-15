#include "fem/CrouzeixRaviartElement.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    CrouzeixRaviartElement::CrouzeixRaviartElement()
      : ScalarElement(FEType::CrouzeixRaviart, 1)
    {
      nr_dof = 3;
      dofs_per_edge = 1;
    }

    // The basis function of edge k is 1 - 2*lambda_j with the vertex j = (k+2)%3
    // opposite to the edge

    double CrouzeixRaviartElement::Value(int k, double x, double y) const
    {
      const double lambda[3] = {1.-x-y, x, y};

      return 1. - 2.*lambda[(k+2)%3];
    }

    Vector2D CrouzeixRaviartElement::Gradient(int k, double, double) const
    {
      const Vector2D grad_lambda[3] = {Vector2D(-1., -1.), Vector2D(1., 0.), Vector2D(0., 1.)};

      return -2.*grad_lambda[(k+2)%3];
    }

    Matrix2D CrouzeixRaviartElement::Hessian(int, double, double) const
    {
      return Matrix2D(0., 0., 0., 0.);
    }

    Coordinate CrouzeixRaviartElement::NodalPoint(int k) const
    {
      double xi, eta;
      EdgeToRefCoords(k, 0.5, xi, eta);

      return Coordinate{xi, eta};
    }

  }
}
