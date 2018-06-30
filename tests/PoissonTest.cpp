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

Vector exact_grad(double x, double y)
{
  Vector grad(2);
  grad[0] = 16.*y*(1.-y)*(1.-2.*x); 
  grad[1] = 16.*x*(1.-x)*(1.-2.*y);
  return grad;
}

int main()
{
  for(int size=10; size<200; size*=2)
    {
      UnitSquareMesh mesh(size);

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

      FEFunction Sol(Space.Interpolate(exact));
      //Sol.CreateFunction(X);
      Sol.WriteVtk("solution.vtk");
    
      ErrorNorm Error;
      Error.SetExactValue(&exact);
      Error.SetExactGradient(&exact_grad);
      Error.SetFEFunction(Sol);
      std::cout << "L2 Error: " << Error.Compute(L2) << std::endl;  
      std::cout << "H1 Error: " << Error.Compute(H1_SEMI) << std::endl;
    }  
  return 0;
}
