#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/DGElement.h"
#include "fem/FEFunction.h"
#include "fem/ErrorNorm.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

// L2 projection onto the DG elements of degree 0 on the unit square, (u_h, v) = (f, v).
// The mass matrix is diagonal here, so u_h is the mean value of f on each cell and the
// error converges with order 1. The space is the pressure space of the lowest order
// Raviart-Thomas discretization.

double smooth(const Coordinate& p)
{
  return exp(p.x)*(1. + p.y*p.y);
}

double constant(const Coordinate&)
{
  return 2.5;
}

/// Projects f onto the space
FEFunction Project(FESpace& Space, ScalarFunction f)
{
  BilinearForm M(Space, Space);
  M.AddVolumeTerm([](const PointValues& u, const PointValues& v)
                  { return u.value * v.value; });
  M.Assemble();

  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();

  FEFunction Sol(Space);
  Sol.CreateFunction(M.SystemMatrix().Solve(F.LoadVector(), LIN_SOLVER::UMFPACK));

  return Sol;
}

int main()
{
  const int levels = 5;

  Mesh mesh = UnitSquareMesh(3);

  std::vector<double> l2_errors;
  bool ok = true;

  // A constant lies in the space, so it has to be reproduced exactly
  {
    DGElement element(0);
    FESpace Space(mesh, element);

    FEFunction Sol = Project(Space, constant);

    ErrorNorm Error(constant);
    const double deviation = Error.Compute(Sol, L2);

    std::cout << "Projection of a constant, L2 deviation " << std::scientific
              << std::setprecision(3) << deviation << std::endl << std::endl;

    if(deviation > 1.e-12)
      {
        std::cerr << "ERROR: a constant is not reproduced exactly.\n";
        ok = false;
      }
  }

  std::cout << std::setw(8) << "Cells" << std::setw(8) << "DOFs"
            << std::setw(12) << "L2 error" << std::setw(8) << "eoc" << std::endl;

  for(int level=0; level<levels; ++level)
    {
      DGElement element(0);
      FESpace Space(mesh, element);

      if(Space.NrDof() != mesh.NrCells())
        {
          std::cerr << "ERROR: " << Space.NrDof() << " DOFs, but the mesh has "
                    << mesh.NrCells() << " cells.\n";
          ok = false;
        }

      FEFunction Sol = Project(Space, smooth);

      ErrorNorm Error(smooth);
      l2_errors.push_back(Error.Compute(Sol, L2));

      std::cout << std::setw(8) << mesh.NrCells() << std::setw(8) << Space.NrDof()
                << std::scientific << std::setprecision(3) << std::setw(12)
                << l2_errors.back();
      if(level > 0)
        std::cout << std::fixed << std::setw(8) << log2(l2_errors[level-1]/l2_errors[level]);
      std::cout << std::endl;

      if(level+1 < levels)
        {
          mesh.RefineUniform();
          mesh.RefineUniform();
        }
    }

  const double l2_eoc = log2(l2_errors[levels-2]/l2_errors[levels-1]);

  if(std::fabs(l2_eoc - 1.) > 0.15)
    {
      std::cerr << "ERROR: the L2 error does not converge with order 1.\n";
      ok = false;
    }

  if(!ok)
    return 1;

  std::cout << "\nDGTest was successful.\n";
  return 0;
}
