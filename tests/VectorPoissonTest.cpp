#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/ProductElement.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "fem/ErrorNorm.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// -Laplace(u) = f  for a vector valued u on the unit square, u = 0 on the boundary. Both
// components are solved in one space [P_k]^2, so the whole problem is a single bilinear
// form with one integrand.

double ux(const Coordinate& p) { return 16.*p.x*(1.-p.x)*p.y*(1.-p.y); }
double uy(const Coordinate& p) { return sin(M_PI*p.x)*sin(M_PI*p.y); }

Vector2D grad_ux(const Coordinate& p)
{
  return Vector2D(16.*(1.-2.*p.x)*p.y*(1.-p.y), 16.*p.x*(1.-p.x)*(1.-2.*p.y));
}

Vector2D grad_uy(const Coordinate& p)
{
  return Vector2D(M_PI*cos(M_PI*p.x)*sin(M_PI*p.y), M_PI*sin(M_PI*p.x)*cos(M_PI*p.y));
}

Vector2D f(const Coordinate& p)
{
  return Vector2D(32.*(p.x*(1.-p.x) + p.y*(1.-p.y)), 2.*M_PI*M_PI*uy(p));
}

/// (grad u, grad v) for vector valued functions
double Stiffness(const VectorValues& u, const VectorValues& v)
{
  return ddot(u.gradient, v.gradient);
}

/// (f, v)
double Force(const QuadPoint& p, const VectorValues& v)
{
  return dot(f(p.x), v.value);
}

bool Run(int degree, int levels)
{
  Mesh mesh = UnitSquareMesh(3);

  std::vector<double> l2_errors, h1_errors;

  std::cout << "\nP" << degree << "\n" << std::string(52, '=') << "\n"
            << std::setw(8) << "Cells" << std::setw(8) << "DOFs"
            << std::setw(12) << "L2 error" << std::setw(8) << "eoc"
            << std::setw(12) << "H1 error" << std::setw(8) << "eoc" << std::endl;

  for(int level=0; level<levels; ++level)
    {
      LagrangeElement scalar(degree);
      ProductElement element(scalar, 2);

      FESpace V(mesh, element);

      BilinearForm A(V, V);
      A.AddVolumeTerm(Stiffness);
      A.Assemble();

      LinearForm F(V);
      F.AddVolumeTerm(Force);
      F.Assemble();

      FEFunction Sol(V);
      Sol.CreateFunction(A.SystemMatrix().Solve(F.LoadVector(), LIN_SOLVER::UMFPACK));

      ErrorNorm Ex(ux, grad_ux), Ey(uy, grad_uy);

      l2_errors.push_back(std::hypot(Ex.Compute(Sol, L2, 0), Ey.Compute(Sol, L2, 1)));
      h1_errors.push_back(std::hypot(Ex.Compute(Sol, H1_SEMI, 0), Ey.Compute(Sol, H1_SEMI, 1)));

      std::cout << std::setw(8) << mesh.NrCells() << std::setw(8) << V.NrFreeDof()
                << std::scientific << std::setprecision(3)
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

  std::cout << "\nVectorPoissonTest was successful.\n";
  return 0;
}
