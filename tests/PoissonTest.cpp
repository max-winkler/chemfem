#include <iostream>
#include <cmath>

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
  return sin(2*x*M_PI)+cos(2*y*M_PI);
}

double c(double x, double y)
{
  return 1.;
}

int main()
{
  UnitSquareMesh mesh(3);

  std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
  std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;
  
  LagrangeElement element(3);
  FESpace Space(mesh, element);

  BilinearForm Laplace(Space, Space);
  Laplace.AddLaplaceTerm();
  Laplace.AddReactionTerm(c);
  Laplace.Assemble();

  SparseMatrix& Matrix = Laplace.SystemMatrix();
  
  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();

  Vector& Vec = F.LoadVector();

  Vector X(Matrix.Solve(Vec));

  Vector Res(Matrix*X - Vec);
  std::cout << "Error of equation system: " << Res.Norm() << std::endl;

  mesh.WriteVtk("solution.vtk", X);
  
  return 0;
}
