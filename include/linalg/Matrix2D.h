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

      double Trace() const { return a00 + a11; }

      Matrix2D Transpose() const { return Matrix2D(a00, a10, a01, a11); }

      Matrix2D& operator+=(const Matrix2D& B)
      { a00 += B.a00; a01 += B.a01; a10 += B.a10; a11 += B.a11; return *this; }

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

    inline Matrix2D operator*(const Matrix2D& A, const Matrix2D& B)
    { return Matrix2D(A.a00*B.a00 + A.a01*B.a10, A.a00*B.a01 + A.a01*B.a11,
                      A.a10*B.a00 + A.a11*B.a10, A.a10*B.a01 + A.a11*B.a11); }

    /// Double contraction A : B, the scalar product of two matrices
    inline double ddot(const Matrix2D& A, const Matrix2D& B)
    { return A.a00*B.a00 + A.a01*B.a01 + A.a10*B.a10 + A.a11*B.a11; }

    inline Matrix2D operator*(double s, const Matrix2D& A)
    { return Matrix2D(s*A.a00, s*A.a01, s*A.a10, s*A.a11); }

  };
};

#endif
