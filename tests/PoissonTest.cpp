#include <iostream>
#include <cmath>

#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

double f(double x, double y)
{
  return 1.;//sin(2*x*M_PI)+cos(2*y*M_PI);
}

double c(double x, double y)
{
  return 1.;
}

int main()
{
  UnitSquareMesh mesh(16);

  std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
  std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;
  
  LagrangeElement element(2);
  FESpace Space(mesh, element);

  BilinearForm Laplace(Space, Space);
  Laplace.AddLaplaceTerm();
  //Laplace.AddReactionTerm(c);
  Laplace.Assemble();

  SparseMatrix& Matrix = Laplace.SystemMatrix();

  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();

  FEFunction Sol(Space);
  
  Vector& Vec = F.LoadVector();
    
  Vector X(Matrix.Solve(Vec));
    
  Vector Res(Matrix*X - Vec);
  std::cout << "Error of equation system: " << Res.Norm() << std::endl;
  
  Vector X_full(Space.IncorporateBC(X));
    
  mesh.WriteVtk("solution.vtk", X_full);
  
  return 0;
}
