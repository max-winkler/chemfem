#ifndef _VECTOR_2D_H_
#define _VECTOR_2D_H_

#include <cmath>

#include "linalg/Coordinate.h"

namespace chemfem{
  namespace linalg{

    /**
     * A vector of the plane, for gradients, normals and displacements.
     */
    struct Vector2D
    {
      double x, y;

      Vector2D() : x(0.), y(0.) {}
      Vector2D(double x, double y) : x(x), y(y) {}

      double operator[](int i) const { return i == 0 ? x : y; }
      double& operator[](int i) { return i == 0 ? x : y; }

      Vector2D& operator+=(const Vector2D& v) { x += v.x; y += v.y; return *this; }
      Vector2D& operator-=(const Vector2D& v) { x -= v.x; y -= v.y; return *this; }
      Vector2D& operator*=(double s) { x *= s; y *= s; return *this; }

      double Norm() const { return std::sqrt(x*x + y*y); }
    };

    inline Vector2D operator+(const Vector2D& a, const Vector2D& b)
    { return Vector2D(a.x + b.x, a.y + b.y); }

    inline Vector2D operator-(const Vector2D& a, const Vector2D& b)
    { return Vector2D(a.x - b.x, a.y - b.y); }

    inline Vector2D operator*(double s, const Vector2D& v)
    { return Vector2D(s*v.x, s*v.y); }

    inline double dot(const Vector2D& a, const Vector2D& b)
    { return a.x*b.x + a.y*b.y; }

    inline Coordinate operator+(const Coordinate& p, const Vector2D& v)
    { return Coordinate{p.x + v.x, p.y + v.y}; }

    inline Coordinate operator-(const Coordinate& p, const Vector2D& v)
    { return Coordinate{p.x - v.x, p.y - v.y}; }

    inline Vector2D operator-(const Coordinate& p, const Coordinate& q)
    { return Vector2D(p.x - q.x, p.y - q.y); }

  };
};

#endif
