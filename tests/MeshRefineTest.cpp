#include <iostream>
#include <cmath>

#include "mesh/UnitSquareMesh.h"
#include "fem/LinearForm.h"
#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "fem/FEFunction.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::mesh;
using namespace chemfem::linalg;
using namespace chemfem::fem;

double f(double x, double y)
{
  return 10.*sin(4*x*M_PI)*cos(2*y*M_PI);
}

double c(double x, double y)
{
  return 1.;
}

int main()
{
  // Build up and refine mesh
  Mesh mesh = UnitSquareMesh(4);
  mesh.WriteVtk("mesh_old.vtk");

  for(int k=0; k<2; ++k)
    {
      mesh = mesh.RefineUniform();  
      mesh.Check();
    }
  mesh.WriteVtk("mesh.vtk");

  // Solve Poisson equation
  LagrangeElement element(1);
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

  FEFunction Sol(Space);
  Sol.CreateFunction(X);
  Sol.WriteVtk("solution.vtk");
  
  return 0;
}
