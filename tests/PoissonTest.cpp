#include <iostream>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

double f(double x, double y)
{
  return 1.;
}

int main()
{
  UnitSquareMesh mesh(4);

  std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
  std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;
  
  LagrangeElement element(1);
  FESpace Space(mesh, element);

  BilinearForm Laplace(Space, Space);
  Laplace.AddLaplaceTerm();
  Laplace.Assemble();

  SparseMatrix& Matrix = Laplace.SystemMatrix();
  
  std::cout << "The finite element system matrix is:\n" << Matrix << std::endl;

  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();

  Vector& Vec = F.LoadVector();

  std::cout << "The load vector is:\n" << Vec << std::endl;

  Vector X(Matrix.Solve(Vec));
  
  return 0;
}
