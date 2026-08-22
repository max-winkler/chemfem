#ifndef _MATRIX_2D_H_
#define _MATRIX_2D_H_

#include "linalg/Vector2D.h"

namespace chemfem{
  namespace linalg{

    /**
     * A 2x2 matrix, for the Jacobian of the reference transformation and its inverse.
     */
    struct Matrix2D
    {
      double a00, a01;
      double a10, a11;

      Matrix2D() : a00(0.), a01(0.), a10(0.), a11(0.) {}
      Matrix2D(double a00, double a01, double a10, double a11)
        : a00(a00), a01(a01), a10(a10), a11(a11) {}

      double Determinant() const { return a00*a11 - a01*a10; }

      Matrix2D Transpose() const { return Matrix2D(a00, a10, a01, a11); }

      Matrix2D Invert() const
      {
        const double s = 1./Determinant();
        return Matrix2D(a11*s, -a01*s, -a10*s, a00*s);
      }
    };

    inline Vector2D operator*(const Matrix2D& A, const Vector2D& v)
    { return Vector2D(A.a00*v.x + A.a01*v.y, A.a10*v.x + A.a11*v.y); }

    inline Vector2D operator*(const Matrix2D& A, const Coordinate& p)
    { return Vector2D(A.a00*p.x + A.a01*p.y, A.a10*p.x + A.a11*p.y); }

  };
};

#endif
