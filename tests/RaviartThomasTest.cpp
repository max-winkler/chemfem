#include <iostream>
#include <iomanip>
#include <cmath>

#include "fem/RaviartThomasElement.h"
#include "fem/PointValues.h"
#include "fem/FEFunction.h"
#include "mesh/UnitSquareMesh.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// The defining property of the lowest order Raviart-Thomas element, checked directly on the
// element without a mesh or an assembly: the flux of the basis function k through the edge j
// is 1 for k == j and 0 otherwise. The second part repeats it after the contravariant Piola
// transform on a concrete triangle, where the fluxes have to come out the same.
//
// The third part is the property the discretization rests on and the one the edge signs are
// there for: on a mesh, the normal component of an FE function is single valued across an
// interior edge, although the tangential one jumps.

const Coordinate Reference[3] = {{0., 0.}, {1., 0.}, {0., 1.}};
const Coordinate Physical[3] = {{0., 0.}, {2., 0.}, {1., 3.}};

bool Nowhere(const Coordinate&)
{
  return false;
}

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

  // Normal continuity on a mesh
  Mesh mesh = UnitSquareMesh(3);
  mesh.RefineUniform();

  RaviartThomasElement element;
  FESpace Space(mesh, element, Nowhere);

  if(Space.NrDof() != mesh.NrEdges())
    {
      std::cerr << "ERROR: " << Space.NrDof() << " DOFs, but the mesh has " << mesh.NrEdges()
                << " edges.\n";
      ok = false;
    }

  // Without a single reversed edge the continuity below would hold for trivial reasons
  size_t negative = 0;
  for(size_t c=0; c<mesh.NrCells(); ++c)
    for(size_t k=0; k<Space.NrLocalDof(); ++k)
      if(Space.LocalSign(c, k) < 0.)
        ++negative;

  std::cout << "Edge signs: " << negative << " of " << mesh.NrCells()*Space.NrLocalDof()
            << " local DOFs carry -1" << std::endl;

  if(negative == 0)
    {
      std::cerr << "ERROR: no local DOF is reversed, the test below proves nothing.\n";
      ok = false;
    }

  // One function that involves every basis function
  Vector Coefficients(Space.NrDof());
  for(size_t d=0; d<Coefficients.size(); ++d)
    Coefficients[d] = sin(1. + d);

  FEFunction Sigma(Space);
  Sigma.SetCoefficients(Coefficients);

  const std::vector<Edge>& Edges = mesh.GetEdgeList();
  const std::vector<Cell>& Cells = mesh.GetCellList();

  double worst = 0.;
  size_t interior = 0;

  for(size_t e=0; e<Edges.size(); ++e)
    {
      if(Edges[e].Type() != EdgeType::INTERFACE_EDGE)
        continue;

      ++interior;

      const std::pair<size_t, size_t> neighbors = Edges[e].GetNeighbors();

      const int a = Cells[neighbors.first].EdgeIndex(e);
      const int b = Cells[neighbors.second].EdgeIndex(e);

      const Vector2D normal = mesh.GetCellInfo(neighbors.first).Normal(a);

      // The midpoint of the edge is the same physical point from both cells, whichever way
      // they run through it
      double xi_a, eta_a, xi_b, eta_b;
      EdgeToRefCoords(a, 0.5, xi_a, eta_a);
      EdgeToRefCoords(b, 0.5, xi_b, eta_b);

      // VectorValue uses the cell and the reference coordinates, not the physical point
      const Vector2D from_a = Sigma.VectorValue(QuadPoint{Coordinate{0., 0.},
                                                          neighbors.first, xi_a, eta_a});
      const Vector2D from_b = Sigma.VectorValue(QuadPoint{Coordinate{0., 0.},
                                                          neighbors.second, xi_b, eta_b});

      worst = std::max(worst, std::fabs(dot(from_a, normal) - dot(from_b, normal)));
    }

  std::cout << "Normal component across " << interior << " interior edges, largest jump "
            << std::scientific << std::setprecision(3) << worst << std::endl;

  if(worst > 1.e-12)
    {
      std::cerr << "ERROR: the normal component is not continuous.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nRaviartThomasTest was successful.\n";
  return 0;
}
