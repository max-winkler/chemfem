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

// -Laplace(u) + b.grad(u) = f  on the unit square, u = 0 on the boundary,
// with the rotational field b = (y,-x) and the exact solution below.

Vector2D b(const Coordinate& p)
{
  return Vector2D(p.y, -p.x);
}

double exact(const Coordinate& p)
{
  return 16.*p.x*(1.-p.x)*p.y*(1.-p.y);
}

Vector2D exact_grad(const Coordinate& p)
{
  return Vector2D(16.*p.y*(1.-p.y)*(1.-2.*p.x),
                  16.*p.x*(1.-p.x)*(1.-2.*p.y));
}

double f(const Coordinate& p)
{
  const double laplace = -32.*(p.x*(1.-p.x) + p.y*(1.-p.y));

  return -laplace + dot(b(p), exact_grad(p));
}

int main()
{
  const int max_iter = 6;

  std::vector<double> l2_errors, h1_errors;
  std::vector<size_t> cells;

  Mesh mesh = UnitSquareMesh(2);

  for(int iter=0; iter<max_iter; ++iter)
    {
      LagrangeElement element(1);
      FESpace Space(mesh, element);

      BilinearForm A(Space, Space);
      A.AddLaplaceTerm();
      A.AddConvectionTerm(b);
      A.Assemble();

      LinearForm F(Space);
      F.AddVolumeForce(f);
      F.Assemble();

      FEFunction Sol(Space);
      Sol.CreateFunction(A.SystemMatrix().Solve(F.LoadVector()));

      ErrorNorm Error;
      Error.SetExactValue(exact);
      Error.SetExactGradient(exact_grad);
      Error.SetFEFunction(Sol);

      cells.push_back(mesh.NrCells());
      l2_errors.push_back(Error.Compute(L2));
      h1_errors.push_back(Error.Compute(H1_SEMI));

      if(iter+1 < max_iter)
        {
          mesh.RefineUniform();
          mesh.RefineUniform();
        }
    }

  std::cout << "\n" << std::setw(10) << "Zellen"
            << std::setw(18) << "L2-error" << std::setw(12) << "L2-eoc"
            << std::setw(18) << "H1-error" << std::setw(12) << "H1-eoc"
            << std::endl;

  bool ok = true;

  for(int iter=0; iter<max_iter; ++iter)
    {
      double l2_eoc = 0., h1_eoc = 0.;

      if(iter > 0)
        {
          l2_eoc = log(l2_errors[iter]/l2_errors[iter-1])/log(0.5);
          h1_eoc = log(h1_errors[iter]/h1_errors[iter-1])/log(0.5);
        }

      std::cout << std::setw(10) << cells[iter]
                << std::setw(18) << l2_errors[iter] << std::setw(12) << l2_eoc
                << std::setw(18) << h1_errors[iter] << std::setw(12) << h1_eoc
                << std::endl;

      // The two finest levels have to show the rates of a P1 discretization
      if(iter >= max_iter-2)
        {
          if(std::fabs(l2_eoc - 2.) > 0.1) ok = false;
          if(std::fabs(h1_eoc - 1.) > 0.1) ok = false;
        }
    }

  if(!ok)
    {
      std::cerr << "\nERROR: the convergence rates are not the expected ones "
                << "(2 in L2 and 1 in H1).\n";
      return 1;
    }

  std::cout << "\nConvectionTest was successful.\n";
  return 0;
}
