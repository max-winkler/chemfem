#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/CrouzeixRaviartElement.h"
#include "fem/FEFunction.h"
#include "fem/ErrorNorm.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// -Laplace(u) = f  on the unit square, u = 0 on the boundary, discretized with the
// nonconforming Crouzeix-Raviart element. The H1 error is the broken one, taken cell
// by cell.

double exact(const Coordinate& p)
{
  return sin(M_PI*p.x)*sin(M_PI*p.y);
}

Vector2D exact_grad(const Coordinate& p)
{
  return Vector2D(M_PI*cos(M_PI*p.x)*sin(M_PI*p.y), M_PI*sin(M_PI*p.x)*cos(M_PI*p.y));
}

double f(const Coordinate& p)
{
  return 2.*M_PI*M_PI*exact(p);
}

int main()
{
  const int levels = 5;

  Mesh mesh = UnitSquareMesh(3);

  std::vector<double> l2_errors, h1_errors;
  bool ok = true;

  std::cout << std::setw(8) << "Cells" << std::setw(8) << "DOFs"
            << std::setw(12) << "L2 error" << std::setw(8) << "eoc"
            << std::setw(12) << "H1 error" << std::setw(8) << "eoc" << std::endl;

  for(int level=0; level<levels; ++level)
    {
      CrouzeixRaviartElement element;
      FESpace Space(mesh, element, WholeBoundary);

      if(Space.NrDof() != mesh.NrEdges())
        {
          std::cerr << "ERROR: " << Space.NrDof() << " DOFs, but the mesh has "
                    << mesh.NrEdges() << " edges.\n";
          ok = false;
        }

      BilinearForm A(Space, Space);
      A.AddLaplaceTerm();
      A.Assemble();

      LinearForm F(Space);
      F.AddVolumeForce(f);
      F.Assemble();

      FEFunction Sol(Space);
      Sol.CreateFunction(A.SystemMatrix().Solve(F.LoadVector(), LIN_SOLVER::UMFPACK));

      ErrorNorm Error;
      Error.SetExactValue(exact);
      Error.SetExactGradient(exact_grad);
      Error.SetFEFunction(Sol);

      l2_errors.push_back(Error.Compute(L2));
      h1_errors.push_back(Error.Compute(H1_SEMI));

      std::cout << std::setw(8) << mesh.NrCells() << std::setw(8) << Space.NrDof()
                << std::scientific << std::setprecision(3) << std::setw(12) << l2_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8) << log2(l2_errors[level-1]/l2_errors[level]);
      else
        std::cout << std::setw(8) << "";
      std::cout << std::scientific << std::setw(12) << h1_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8) << log2(h1_errors[level-1]/h1_errors[level]);
      std::cout << std::endl;

      if(level+1 < levels)
        {
          mesh.RefineUniform();
          mesh.RefineUniform();
        }
    }

  const double l2_eoc = log2(l2_errors[levels-2]/l2_errors[levels-1]);
  const double h1_eoc = log2(h1_errors[levels-2]/h1_errors[levels-1]);

  if(std::fabs(l2_eoc - 2.) > 0.15)
    {
      std::cerr << "ERROR: the L2 error does not converge with order 2.\n";
      ok = false;
    }

  if(std::fabs(h1_eoc - 1.) > 0.15)
    {
      std::cerr << "ERROR: the broken H1 error does not converge with order 1.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nCrouzeixRaviartTest was successful.\n";
  return 0;
}
