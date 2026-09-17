#include <iostream>
#include <iomanip>
#include <cmath>

#include "fem/LinearForm.h"
#include "fem/LinearSystem.h"
#include "fem/BilinearForm.h"
#include "fem/DirichletValues.h"
#include "fem/LagrangeElement.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// A system may be assembled more than once. A Newton iteration does exactly that: it
// rebuilds the matrix from the current iterate in every step. Assembling twice without
// changing anything in between has to give the same right hand side, so the lifting of the
// prescribed values must not accumulate.

double Prescribed(const Coordinate& p)
{
  return p.x + 2.*p.y;
}

double f(const Coordinate&)
{
  return 1.;
}

int main()
{
  Mesh mesh = UnitSquareMesh(3);
  mesh.RefineUniform();

  LagrangeElement P1(1);
  FESpace V(mesh, P1, WholeBoundary);

  DirichletValues g(V);
  g.Set(Prescribed);

  BilinearForm A(V, V);
  A.AddLaplaceTerm();

  LinearForm F(V);
  F.AddVolumeForce(f);

  LinearSystem S(V);
  S.AddLhs(A);
  S.AddRhs(F);
  S.SetDirichletValues(g);

  S.AssembleMatrix();
  const Vector FirstRhs = S.AssembleRhs();
  const Vector FirstSolution = S.Solve(FirstRhs);

  S.AssembleMatrix();
  const Vector SecondRhs = S.AssembleRhs();
  const Vector SecondSolution = S.Solve(SecondRhs);

  const double rhs_error = (FirstRhs - SecondRhs).Norm() / FirstRhs.Norm();
  const double solution_error
    = (FirstSolution - SecondSolution).Norm() / FirstSolution.Norm();

  std::cout << std::scientific << std::setprecision(3)
            << "Right hand side of the second assembly, relative deviation " << rhs_error
            << "\nSolution of the second assembly, relative deviation      "
            << solution_error << std::endl;

  bool ok = true;

  if(rhs_error > 1.e-12)
    {
      std::cerr << "ERROR: assembling twice changes the right hand side. The lifting of "
                << "the prescribed values accumulates.\n";
      ok = false;
    }

  if(solution_error > 1.e-12)
    {
      std::cerr << "ERROR: assembling twice changes the solution.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nReassembleTest was successful.\n";
  return 0;
}
