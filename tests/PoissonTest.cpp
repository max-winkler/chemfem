#include <iostream>
#include <iomanip>
#include <fstream>
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

double f(const Coordinate& p)
{
  return 32.*(p.x*(1.-p.x) + p.y*(1.-p.y)) + 16.*p.x*(1.-p.x)*p.y*(1.-p.y);
}

double c(const Coordinate& p)
{
  return 1.;
}

double exact(const Coordinate& p)
{
  return 16.*p.x*(1.-p.x)*p.y*(1.-p.y);
}

Vector2D exact_grad(const Coordinate& p)
{
  Vector2D grad;
  grad[0] = 16.*p.y*(1.-p.y)*(1.-2.*p.x); 
  grad[1] = 16.*p.x*(1.-p.x)*(1.-2.*p.y);
  return grad;
}

int main()
{
  Mesh mesh = UnitSquareMesh(2);
  
  const int max_iter = 8;

  std::vector<double> l2_errors, h1_errors;
  
  for(int iter=0; iter<max_iter; ++iter)
    {
      std::cout << " Computation on level " << iter+1 << std::endl;
      std::cout << "======================================\n";
      
      std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
      std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;

      // std::cout << mesh << std::endl;
      
      LagrangeElement element(1);
      FESpace Space(mesh, element);
      
      std::cout << "Space:\n" << Space << std::endl;
      
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
      
      FEFunction Sol(Space);
      Sol.CreateFunction(X);
      Sol.WriteVtk("solution.vtk");
    
      ErrorNorm Error;
      Error.SetExactValue(&exact);
      Error.SetExactGradient(&exact_grad);
      Error.SetFEFunction(Sol);

      double l2_error = Error.Compute(L2);
      double h1_error = Error.Compute(H1_SEMI);
      
      l2_errors.push_back(l2_error);
      h1_errors.push_back(h1_error);
      
      if(iter+1 < max_iter)
        {
          mesh.RefineUniform();
          mesh.RefineUniform();
        }
    }  

  std::cout << std::setw(10) << "Iteration"
	  << std::setw(20) << "L2-error" << std::setw(20) << "L2-eoc"
	  << std::setw(20) << "H1-error" << std::setw(20) << "H1-eoc"
	  << std::endl;
  
  for(int iter=0; iter<max_iter; ++iter)
    {
      double l2_eoc = 0., h1_eoc = 0.;
      
      if(iter>0)
	{
	  l2_eoc = log(l2_errors[iter] / l2_errors[iter-1]) / log(0.5);
	  h1_eoc = log(h1_errors[iter] / h1_errors[iter-1]) / log(0.5);
	}

      std::cout << std::setw(10) << iter
		<< std::setw(20) << l2_errors[iter] << std::setw(20) << l2_eoc
		<< std::setw(20) << h1_errors[iter] << std::setw(20) << h1_eoc
		<< std::endl;
    }
  return 0;
}
