#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// ---------------------------------------------------------------------------------
// The three solvers behind SparseMatrix::Solve have to agree on the same system.
//
// The Poisson matrix is symmetric and positive definite, so all three apply. The
// convection term makes it non symmetric, which is where CG loses its justification
// and where a wrong transpose in the UMFPACK binding would show.
// ---------------------------------------------------------------------------------

Vector2D wind(const Coordinate& p)
{
  return Vector2D(p.y, -p.x);
}

double f(const Coordinate& p)
{
  return 1. + p.x + p.y;
}

double Norm2(const Vector& v)
{
  return sqrt(dot(v, v));
}

/// ||A x - b|| / ||b||
double Residual(SparseMatrix& A, const Vector& x, const Vector& b)
{
  return Norm2(A*x - b) / Norm2(b);
}

/// ||x - y|| / ||x||
double Deviation(const Vector& x, const Vector& y)
{
  return Norm2(x - y) / Norm2(x);
}

bool CheckSystem(const std::string& caption, bool convection, bool with_cg)
{
  Mesh mesh = UnitSquareMesh(4);

  LagrangeElement element(2);
  FESpace Space(mesh, element);

  BilinearForm A(Space, Space);
  A.AddLaplaceTerm();
  if(convection)
    A.AddConvectionTerm(wind);
  A.Assemble();

  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();

  SparseMatrix& Matrix = A.SystemMatrix();
  Vector& Rhs = F.LoadVector();

  Vector x_umf(Matrix.Solve(Rhs, LIN_SOLVER::UMFPACK));
  Vector x_gmres(Matrix.Solve(Rhs, LIN_SOLVER::GMRES));

  std::cout << "\n" << caption << "\n" << std::string(60, '=') << "\n"
            << std::setw(10) << "Solver" << std::setw(16) << "rel. residual"
            << std::setw(20) << "dev. from UMFPACK" << std::endl;

  std::cout << std::scientific << std::setprecision(3)
            << std::setw(10) << "UMFPACK" << std::setw(16) << Residual(Matrix, x_umf, Rhs)
            << std::setw(20) << "-" << std::endl
            << std::setw(10) << "GMRES" << std::setw(16) << Residual(Matrix, x_gmres, Rhs)
            << std::setw(20) << Deviation(x_umf, x_gmres) << std::endl;

  bool ok = true;

  // The direct solver is only limited by round off
  if(Residual(Matrix, x_umf, Rhs) > 1.e-10)
    {
      std::cerr << "ERROR: UMFPACK left a residual.\n";
      ok = false;
    }

  if(Deviation(x_umf, x_gmres) > 1.e-6)
    {
      std::cerr << "ERROR: GMRES and UMFPACK disagree.\n";
      ok = false;
    }

  if(with_cg)
    {
      Vector x_cg(Matrix.Solve(Rhs, LIN_SOLVER::CG));

      std::cout << std::setw(10) << "CG" << std::setw(16) << Residual(Matrix, x_cg, Rhs)
                << std::setw(20) << Deviation(x_umf, x_cg) << std::endl;

      if(Deviation(x_umf, x_cg) > 1.e-6)
        {
          std::cerr << "ERROR: CG and UMFPACK disagree.\n";
          ok = false;
        }
    }

  std::cout << std::defaultfloat << std::setprecision(6);

  return ok;
}

int main()
{
  bool ok = true;

  // CG only applies to the symmetric system
  ok &= CheckSystem("Poisson, symmetric and positive definite", false, true);
  ok &= CheckSystem("Convection diffusion, non symmetric", true, false);

  if(!ok)
    return 1;

  std::cout << "\nSolverTest was successful.\n";
  return 0;
}
