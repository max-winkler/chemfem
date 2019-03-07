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

double f(double x, double y)
{
  return 32.*(x*(1.-x) + y*(1.-y)) + 16.*x*(1.-x)*y*(1.-y);
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
  Mesh mesh = UnitSquareMesh(5);
  const int max_iter = 3;

  std::vector<double> l2_errors, h1_errors;
  
  for(int iter=0; iter<max_iter; ++iter)
    {
      std::cout << " Computation on level " << iter << std::endl;
      std::cout << "======================================\n";
      
      std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
      std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;

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
	mesh = mesh.RefineUniform();
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
