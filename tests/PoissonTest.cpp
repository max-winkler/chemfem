#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "fem/LinearForm.h"
#include "fem/LinearSystem.h"
#include "fem/BilinearForm.h"
#include "fem/DirichletValues.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "fem/ErrorNorm.h"
#include "fem/VtkOutput.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// -Laplace(u) + u = f  on the unit square with the exact solution u = x(1-x)(1-y). It
// vanishes on three sides and equals x(1-x) on the bottom one, so the test also covers
// inhomogeneous Dirichlet values.

double exact(const Coordinate& p)
{
  return p.x*(1.-p.x)*(1.-p.y);
}

Vector2D exact_grad(const Coordinate& p)
{
  return Vector2D((1.-2.*p.x)*(1.-p.y), -p.x*(1.-p.x));
}

double f(const Coordinate& p)
{
  return 2.*(1.-p.y) + exact(p);
}

double One(const Coordinate&)
{
  return 1.;
}

/// The bottom side y = 0, where the solution does not vanish
bool Bottom(const Coordinate& p)
{
  return p.y < 1.e-12;
}

int main()
{
  const int levels = 6;

  Mesh mesh = UnitSquareMesh(2);

  std::vector<double> l2_errors, h1_errors;

  std::cout << std::setw(8) << "Cells" << std::setw(8) << "DOFs"
            << std::setw(12) << "L2 error" << std::setw(8) << "eoc"
            << std::setw(12) << "H1 error" << std::setw(8) << "eoc" << std::endl;

  for(int level=0; level<levels; ++level)
    {
      LagrangeElement element(1);
      FESpace Space(mesh, element, WholeBoundary);

      DirichletValues g(Space);
      g.Set(exact, Bottom);

      BilinearForm A(Space, Space);
      A.AddLaplaceTerm();
      A.AddReactionTerm(One);

      LinearForm F(Space);
      F.AddVolumeForce(f);

      LinearSystem S(Space);
      S.AddLhs(A);
      S.AddRhs(F);
      S.SetDirichletValues(g);
      S.AssembleMatrix();

      FEFunction Sol = S.Extract(S.Solve(S.AssembleRhs()));

      ErrorNorm Error(exact, exact_grad);

      l2_errors.push_back(Error.Compute(Sol, L2));
      h1_errors.push_back(Error.Compute(Sol, H1_SEMI));

      std::cout << std::setw(8) << mesh.NrCells() << std::setw(8) << Space.NrFreeDof()
                << std::scientific << std::setprecision(3)
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

      if(level+1 == levels)
        {
          FEFunction Interpolant = Space.Interpolate(exact);

          VtkOutput out(mesh);
          out.AddScalar("u", Sol);
          out.AddScalar("u_exact", Interpolant);
          out.Write("poisson.vtk");
        }

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

  if(std::fabs(l2_eoc - 2.) > 0.15)
    {
      std::cerr << "ERROR: the L2 error does not converge with order 2.\n";
      ok = false;
    }

  if(std::fabs(h1_eoc - 1.) > 0.15)
    {
      std::cerr << "ERROR: the H1 error does not converge with order 1.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nPoissonTest was successful.\n";
  return 0;
}
