#include <iostream>
#include <cmath>

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

double f(double x, double y)
{
  return 32.*(x*(1.-x) + y*(1.-y));
}

double c(double x, double y)
{
  return 1.;
}

double exact(double x, double y)
{
  return 16.*x*(1.-x)*y*(1.-y);
}

int main()
{
  UnitSquareMesh mesh(50);

  std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
  std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;
  
  LagrangeElement element(1);
  FESpace Space(mesh, element);

  BilinearForm Laplace(Space, Space);
  Laplace.AddLaplaceTerm();
  //Laplace.AddReactionTerm(c);
  Laplace.Assemble();

  SparseMatrix& Matrix = Laplace.SystemMatrix();
  
  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();
  
  Vector& Vec = F.LoadVector();
    
  Vector X(Matrix.Solve(Vec));
    
  Vector Res(Matrix*X - Vec);
  std::cout << "Error of equation system: " << Res.Norm() << std::endl;

  FEFunction Sol(Space);//.Interpolate(exact);
  Sol.CreateFunction(X);
  Sol.WriteVtk("solution.vtk");
    
  ErrorNorm Error;
  Error.SetExactValue(&exact);
  Error.SetFEFunction(Sol);
  std::cout << "Error: " << Error.Compute(L2) << std::endl;  

  return 0;
}
