#include "fem/HermiteElement.h"

#include <cmath>
#include <iostream>

#include "fem/Element.h"
#include "linalg/Matrix2D.h"
#include "linalg/Vector2D.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::Vector2D;

namespace chemfem {
namespace fem {

HermiteElement::HermiteElement() : ScalarElement(FEType::Hermite, 3) {
  dofs_per_vertex = 3;
  dofs_interior = 1;
  nr_dof = 10;
}

Coordinate HermiteElement::NodalPoint(int i) const {
  const Coordinate Vertex[3] = {{0., 0.}, {1., 0.}, {0., 1.}};

  if (i < 9) return Vertex[i / 3];
  if (i == 9) return Coordinate{1. / 3, 1. / 3};

  std::cerr << "Requested nodal point " << i << ", but the element has only "
            << nr_dof << " degrees of freedom.\n";

  return Coordinate{0., 0.};
}

DofDescriptor HermiteElement::Dof(int k) const {
  if (k == 9) return DofDescriptor{DofType::PointValue, -1, -1, 0};

  const int v = k / 3, j = k % 3;

  if (j == 0) return DofDescriptor{DofType::PointValue, v, -1, 0};

  // j == 1 is the edge from v to v+1, j == 2 the one from v to v+2. Both are named by
  // the local edge they lie on, which EdgeVector orients away from v.
  return DofDescriptor{DofType::EdgeDirectionalDerivative, v,
                       j == 1 ? v : (v + 2) % 3, j};
}

double HermiteElement::Value(int i, double x, double y) const {
  switch (i) {
    case 0:
      return 2 * pow(x, 3) + 13 * pow(x, 2) * y - 3 * pow(x, 2) +
             13 * x * pow(y, 2) - 13 * x * y + 2 * pow(y, 3) - 3 * pow(y, 2) +
             1;
    case 1:
      return pow(x, 3) + 3 * pow(x, 2) * y - 2 * pow(x, 2) + 2 * x * pow(y, 2) -
             3 * x * y + x;
    case 2:
      return 2 * pow(x, 2) * y + 3 * x * pow(y, 2) - 3 * x * y + pow(y, 3) -
             2 * pow(y, 2) + y;
    case 3:
      return -2 * pow(x, 3) + 7 * pow(x, 2) * y + 3 * pow(x, 2) +
             7 * x * pow(y, 2) - 7 * x * y;
    case 4:
      return 2 * pow(x, 2) * y + x * pow(y, 2) - x * y;
    case 5:
      return -pow(x, 3) + pow(x, 2) + x * pow(y, 2) - x * y;
    case 6:
      return 7 * pow(x, 2) * y + 7 * x * pow(y, 2) - 7 * x * y - 2 * pow(y, 3) +
             3 * pow(y, 2);
    case 7:
      return pow(x, 2) * y - x * y - pow(y, 3) + pow(y, 2);
    case 8:
      return pow(x, 2) * y + 2 * x * pow(y, 2) - x * y;
    case 9:
      return -27 * pow(x, 2) * y - 27 * x * pow(y, 2) + 27 * x * y;
    default:
      std::cerr << "Requested function value of trial function " << i
		<< ", but the element has only " << nr_dof << " degrees of freedom.\n";
      return 0;
  }
}

Vector2D HermiteElement::Gradient(int i, double x, double y) const {
  if(x > 1 || y > 1 || x < 0 || y < 0 || x+y > 1)
    std::cerr << "Invalid quadrature point.\n";
  switch (i) {
    case 0:
      return Vector2D(
          6 * pow(x, 2) + 26 * x * y - 6 * x + 13 * pow(y, 2) - 13 * y,
          13 * pow(x, 2) + 26 * x * y - 13 * x + 6 * pow(y, 2) - 6 * y);
    case 1:
      return Vector2D(
          3 * pow(x, 2) + 6 * x * y - 4 * x + 2 * pow(y, 2) - 3 * y + 1,
          3 * pow(x, 2) + 4 * x * y - 3 * x);
    case 2:
      return Vector2D(
          4 * x * y + 3 * pow(y, 2) - 3 * y,
          2 * pow(x, 2) + 6 * x * y - 3 * x + 3 * pow(y, 2) - 4 * y + 1);
    case 3:
      return Vector2D(
          -6 * pow(x, 2) + 14 * x * y + 6 * x + 7 * pow(y, 2) - 7 * y,
          7 * pow(x, 2) + 14 * x * y - 7 * x);
    case 4:
      return Vector2D(4 * x * y + pow(y, 2) - y, 2 * pow(x, 2) + 2 * x * y - x);
    case 5:
      return Vector2D(-3 * pow(x, 2) + 2 * x + pow(y, 2) - y, 2 * x * y - x);
    case 6:
      return Vector2D(
          14 * x * y + 7 * pow(y, 2) - 7 * y,
          7 * pow(x, 2) + 14 * x * y - 7 * x - 6 * pow(y, 2) + 6 * y);
    case 7:
      return Vector2D(2 * x * y - y, pow(x, 2) - x - 3 * pow(y, 2) + 2 * y);
    case 8:
      return Vector2D(2 * x * y + 2 * pow(y, 2) - y, pow(x, 2) + 4 * x * y - x);
    case 9:
      return Vector2D(-54 * x * y - 27 * pow(y, 2) + 27 * y,
                      -27 * pow(x, 2) - 54 * x * y + 27 * x);
    default:
      return Vector2D(0., 0.);
  }
}

Matrix2D HermiteElement::Hessian(int i, double x, double y) const {
  switch (i) {
    case 0:
      return Matrix2D(2 * (6 * x + 13 * y - 3), 13 * (2 * x + 2 * y - 1),
                      13 * (2 * x + 2 * y - 1), 2 * (13 * x + 6 * y - 3));
    case 1:
      return Matrix2D(2 * (3 * x + 3 * y - 2), 6 * x + 4 * y - 3,
                      6 * x + 4 * y - 3, 4 * x);
    case 2:
      return Matrix2D(4 * y, 4 * x + 6 * y - 3, 4 * x + 6 * y - 3,
                      2 * (3 * x + 3 * y - 2));
    case 3:
      return Matrix2D(2 * (-6 * x + 7 * y + 3), 7 * (2 * x + 2 * y - 1),
                      7 * (2 * x + 2 * y - 1), 14 * x);
    case 4:
      return Matrix2D(4 * y, 4 * x + 2 * y - 1, 4 * x + 2 * y - 1, 2 * x);
    case 5:
      return Matrix2D(2 * (1 - 3 * x), 2 * y - 1, 2 * y - 1, 2 * x);
    case 6:
      return Matrix2D(14 * y, 7 * (2 * x + 2 * y - 1), 7 * (2 * x + 2 * y - 1),
                      2 * (7 * x - 6 * y + 3));
    case 7:
      return Matrix2D(2 * y, 2 * x - 1, 2 * x - 1, 2 * (1 - 3 * y));
    case 8:
      return Matrix2D(2 * y, 2 * x + 4 * y - 1, 2 * x + 4 * y - 1, 4 * x);
    case 9:
      return Matrix2D(-54 * y, 27 * (-2 * x - 2 * y + 1),
                      27 * (-2 * x - 2 * y + 1), -54 * x);
    default:
      return Matrix2D(0, 0, 0, 0);
  }
}


}
}
