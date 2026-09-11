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

// -Laplace(u) = f  on the unit square with u = sin(x) cos(y). The solution vanishes on
// the left side x = 0, which is the Dirichlet boundary. The other three sides carry
// the Neumann condition grad(u).n = g.

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

bool dirichlet(const Coordinate& p)
{
  return p.x < 1.e-12;
}

// The quadrature points lie inside the edges, so g is never evaluated in a corner
double g(const Coordinate& p)
{
  const Vector2D grad = exact_grad(p);

  if(p.x > 1.-1.e-12)
    return grad.x;
  if(p.y < 1.e-12)
    return -grad.y;
  return grad.y;
}

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
      FESpace Space(mesh, element, dirichlet);

      BilinearForm A(Space, Space);
      A.AddLaplaceTerm();
      A.Assemble();

      LinearForm F(Space);
      F.AddVolumeForce(f);
      F.AddNeumannBC(g);
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
        std::cout << std::fixed << std::setw(8)
                  << log2(l2_errors[level-1]/l2_errors[level]);
      else
        std::cout << std::setw(8) << "";
      std::cout << std::scientific << std::setw(12) << h1_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8)
                  << log2(h1_errors[level-1]/h1_errors[level]);
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

  std::cout << "\nNeumannTest was successful.\n";
  return 0;
}
