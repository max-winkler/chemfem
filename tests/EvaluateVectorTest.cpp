#include <iostream>
#include <iomanip>
#include <cmath>

#include "fem/ProductElement.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "mesh/UnitSquareMesh.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// EvaluateVector returns value, gradient, divergence and rotation of a vector valued FE
// function. A quadratic field lies in [P2]^2 without an error, so the interpolant is the
// field itself and the derivatives have to come out exactly.

Vector2D Field(const Coordinate& p)
{
  return Vector2D(p.x*p.x + 2.*p.x*p.y, 3.*p.y*p.y - p.x);
}

/// Row i is the gradient of the component i
Matrix2D FieldGradient(const Coordinate& p)
{
  return Matrix2D(2.*p.x + 2.*p.y, 2.*p.x,
                  -1.,             6.*p.y);
}

int main()
{
  Mesh mesh = UnitSquareMesh(3);
  mesh.RefineUniform();

  LagrangeElement P2(2);
  ProductElement Velocity(P2, 2);

  FESpace V(mesh, Velocity);

  FEFunction U = V.Interpolate(Field);

  double value_error = 0., gradient_error = 0.;
  double divergence_error = 0., curl_error = 0., cache_error = 0.;

  for(size_t c=0; c<mesh.NrCells(); ++c)
    {
      // The barycenter of the cell is the point xi = eta = 1/3 of its reference element
      const Coordinate x = mesh.GetCellInfo(c).Barycenter();
      const QuadPoint p{x, c, 1./3., 1./3.};

      const VectorValues u = U.EvaluateVector(p);

      const Vector2D value = Field(x);
      const Matrix2D gradient = FieldGradient(x);

      value_error = std::max(value_error, (u.value - value).Norm());

      gradient_error = std::max(gradient_error,
                                std::fabs(u.gradient.a00 - gradient.a00));
      gradient_error = std::max(gradient_error,
                                std::fabs(u.gradient.a01 - gradient.a01));
      gradient_error = std::max(gradient_error,
                                std::fabs(u.gradient.a10 - gradient.a10));
      gradient_error = std::max(gradient_error,
                                std::fabs(u.gradient.a11 - gradient.a11));

      divergence_error = std::max(divergence_error,
                                  std::fabs(u.divergence - gradient.Trace()));

      // The scalar rotation is du_1/dx - du_0/dy
      curl_error = std::max(curl_error,
                            std::fabs(u.curl - (gradient.a10 - gradient.a01)));

      // VectorValue is served from the cache EvaluateVector has just filled
      cache_error = std::max(cache_error, (U.VectorValue(p) - u.value).Norm());
    }

  std::cout << std::scientific << std::setprecision(3)
            << "Value       " << value_error
            << "\nGradient    " << gradient_error
            << "\nDivergence  " << divergence_error
            << "\nRotation    " << curl_error
            << "\nCache       " << cache_error << std::endl;

  bool ok = true;

  const double tolerance = 1.e-12;

  if(value_error > tolerance)
    {
      std::cerr << "ERROR: the value of the interpolant is wrong.\n";
      ok = false;
    }

  if(gradient_error > tolerance)
    {
      std::cerr << "ERROR: the gradient is wrong. Row i has to hold the gradient of the "
                << "component i.\n";
      ok = false;
    }

  if(divergence_error > tolerance)
    {
      std::cerr << "ERROR: the divergence is wrong.\n";
      ok = false;
    }

  if(curl_error > tolerance)
    {
      std::cerr << "ERROR: the rotation is wrong.\n";
      ok = false;
    }

  if(cache_error > tolerance)
    {
      std::cerr << "ERROR: VectorValue and EvaluateVector disagree in the same point.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nEvaluateVectorTest was successful.\n";
  return 0;
}
