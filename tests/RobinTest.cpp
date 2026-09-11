#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "fem/ErrorNorm.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// -Laplace(u) = f  on the unit square with u = sin(x) cos(y). The left side x = 0 is a
// Dirichlet boundary, where u vanishes. The other three sides carry the Robin condition
//
//   grad(u).n + alpha u = g
//
// with a variable coefficient alpha. The weak form has a boundary integral on both sides,
//
//   (grad u, grad v) + (alpha u, v)_Robin = (f, v) + (g, v)_Robin

double exact(const Coordinate& p)
{
  return sin(p.x)*cos(p.y);
}

Vector2D exact_grad(const Coordinate& p)
{
  return Vector2D(cos(p.x)*cos(p.y), -sin(p.x)*sin(p.y));
}

double f(const Coordinate& p)
{
  return 2.*sin(p.x)*cos(p.y);
}

double alpha(const Coordinate& p)
{
  return 1. + p.y;
}

bool Dirichlet(const Coordinate& p)
{
  return p.x < 1.e-12;
}

bool Robin(const Coordinate& p)
{
  return !Dirichlet(p);
}

/// (alpha u, v) on the Robin boundary
struct RobinTerm
{
  double operator()(const QuadPoint& p, const CellGeometry&, const EdgeGeometry&,
                    const PointValues& u, const PointValues& v) const
  {
    return alpha(p.x) * u.value * v.value;
  }
};

/// (g, v) on the Robin boundary, with g = grad(u).n + alpha u and n the normal of the edge
struct RobinData
{
  double operator()(const QuadPoint& p, const CellGeometry&, const EdgeGeometry& edge,
                    const PointValues& v) const
  {
    const double g = dot(exact_grad(p.x), edge.normal) + alpha(p.x)*exact(p.x);
    return g * v.value;
  }
};

bool Run(int degree, int levels)
{
  Mesh mesh = UnitSquareMesh(3);

  std::vector<double> l2_errors, h1_errors;

  std::cout << "\nP" << degree << "\n" << std::string(52, '=') << "\n"
            << std::setw(8) << "Cells" << std::setw(12) << "L2 error" << std::setw(8) << "eoc"
            << std::setw(12) << "H1 error" << std::setw(8) << "eoc" << std::endl;

  for(int level=0; level<levels; ++level)
    {
      LagrangeElement element(degree);
      FESpace Space(mesh, element, Dirichlet);

      BilinearForm A(Space, Space);
      A.AddLaplaceTerm();
      A.AddBoundaryTerm(RobinTerm(), Robin);
      A.Assemble();

      LinearForm F(Space);
      F.AddVolumeForce(f);
      F.AddBoundaryTerm(RobinData(), Robin);
      F.Assemble();

      FEFunction Sol(Space);
      Sol.CreateFunction(A.SystemMatrix().Solve(F.LoadVector(), LIN_SOLVER::UMFPACK));

      ErrorNorm Error;
      Error.SetExactValue(exact);
      Error.SetExactGradient(exact_grad);
      Error.SetFEFunction(Sol);

      l2_errors.push_back(Error.Compute(L2));
      h1_errors.push_back(Error.Compute(H1_SEMI));

      std::cout << std::setw(8) << mesh.NrCells() << std::scientific << std::setprecision(3)
                << std::setw(12) << l2_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8) << log2(l2_errors[level-1]/l2_errors[level]);
      else
        std::cout << std::setw(8) << "";
      std::cout << std::scientific << std::setw(12) << h1_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8) << log2(h1_errors[level-1]/h1_errors[level]);
      std::cout << std::endl;

      // Each call bisects every cell once, so two of them halve h
      if(level+1 < levels)
        {
          mesh.RefineUniform();
          mesh.RefineUniform();
        }
    }

  std::cout << std::defaultfloat << std::setprecision(6);

  const double l2_eoc = log2(l2_errors[levels-2]/l2_errors[levels-1]);
  const double h1_eoc = log2(h1_errors[levels-2]/h1_errors[levels-1]);

  bool ok = true;

  if(std::fabs(l2_eoc - (degree+1)) > 0.15)
    {
      std::cerr << "ERROR: the L2 error does not converge with order " << degree+1 << ".\n";
      ok = false;
    }

  if(std::fabs(h1_eoc - degree) > 0.15)
    {
      std::cerr << "ERROR: the H1 error does not converge with order " << degree << ".\n";
      ok = false;
    }

  return ok;
}

int main()
{
  bool ok = true;

  ok &= Run(1, 5);
  ok &= Run(2, 4);

  if(!ok)
    return 1;

  std::cout << "\nRobinTest was successful.\n";
  return 0;
}
