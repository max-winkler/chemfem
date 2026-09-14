#include <iostream>
#include <iomanip>
#include <cmath>

#include "fem/RaviartThomasElement.h"
#include "fem/PointValues.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;

// The defining property of the lowest order Raviart-Thomas element, checked directly on the
// element without a mesh or an assembly: the flux of the basis function k through the edge j
// is 1 for k == j and 0 otherwise. The second part repeats it after the contravariant Piola
// transform on a concrete triangle, where the fluxes have to come out the same. That is the
// property the whole discretization rests on, because it makes the normal component single
// valued between two cells.

const Coordinate Reference[3] = {{0., 0.}, {1., 0.}, {0., 1.}};
const Coordinate Physical[3] = {{0., 0.}, {2., 0.}, {1., 3.}};

/// Outward unit normal and length of the edge from V[j] to V[j+1] of a positive triangle
void OutwardNormal(const Coordinate* V, int j, Vector2D& normal, double& length)
{
  const Vector2D t(V[(j+1)%3].x - V[j].x, V[(j+1)%3].y - V[j].y);

  length = hypot(t.x, t.y);
  normal = Vector2D(t.y/length, -t.x/length);
}

int main()
{
  const RaviartThomasElement E;

  // Two point Gauss formula on [0,1], exact up to cubic polynomials
  const double offset = 0.5/sqrt(3.);
  const double s[2] = {0.5 - offset, 0.5 + offset};

  const Matrix2D Jac(Physical[1].x - Physical[0].x, Physical[2].x - Physical[0].x,
                     Physical[1].y - Physical[0].y, Physical[2].y - Physical[0].y);

  const double det = (Physical[1].x - Physical[0].x)*(Physical[2].y - Physical[0].y)
    - (Physical[2].x - Physical[0].x)*(Physical[1].y - Physical[0].y);

  bool ok = true;

  std::cout << "Determinant of the reference map " << det << std::endl << std::endl;

  for(int mapped=0; mapped<2; ++mapped)
    {
      const Coordinate* V = mapped ? Physical : Reference;

      std::cout << (mapped ? "Fluxes after the Piola transform" : "Fluxes on the reference"
                    " element") << std::endl;
      std::cout << std::setw(14) << "" << std::setw(12) << "edge 0" << std::setw(12)
                << "edge 1" << std::setw(12) << "edge 2" << std::endl;

      for(int k=0; k<3; ++k)
        {
          std::cout << std::setw(14) << ("basis " + std::to_string(k));

          for(int j=0; j<3; ++j)
            {
              Vector2D normal;
              double length;
              OutwardNormal(V, j, normal, length);

              double flux = 0.;

              for(int q=0; q<2; ++q)
                {
                  double xi, eta;
                  EdgeToRefCoords(j, s[q], xi, eta);

                  const VectorRefValues ref = E.VectorReference(k, xi, eta);

                  const Vector2D value = mapped
                    ? MapFromReference(ref, Jac, det).value : ref.value;

                  flux += 0.5 * dot(value, normal) * length;
                }

              std::cout << std::scientific << std::setprecision(3) << std::setw(12) << flux;

              if(std::fabs(flux - (k == j ? 1. : 0.)) > 1.e-13)
                {
                  std::cerr << "\nERROR: the flux of basis function " << k << " through edge "
                            << j << " is " << flux << ".\n";
                  ok = false;
                }
            }

          std::cout << std::endl;
        }

      std::cout << std::endl;
    }

  // The Piola transform divides the divergence by the determinant
  for(int k=0; k<3; ++k)
    {
      const VectorRefValues ref = E.VectorReference(k, 0.25, 0.25);
      const double divergence = MapFromReference(ref, Jac, det).divergence;

      if(std::fabs(ref.divergence - 2.) > 1.e-14
         || std::fabs(divergence - 2./det) > 1.e-14)
        {
          std::cerr << "ERROR: basis function " << k << " has reference divergence "
                    << ref.divergence << " and mapped divergence " << divergence << ".\n";
          ok = false;
        }
    }

  if(!ok)
    return 1;

  std::cout << "RaviartThomasTest was successful.\n";
  return 0;
}
