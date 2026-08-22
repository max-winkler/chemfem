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

double f(const Coordinate& p)
{
  return 10.*sin(2*p.x*M_PI)*cos(2*p.y*M_PI);
}

double c(const Coordinate& p)
{
  return 1.;
}

int main()
{
  // Build up and refine mesh
  UnitSquareMesh mesh(3);
  mesh.WriteVtk("mesh_old.vtk");

  for(int k=0; k<5; ++k)
    {
      std::cout << "Refinement #" << k << std::endl;
      mesh.RefineUniform();  
      mesh.Check();
    }
  mesh.WriteVtk("mesh.vtk");

  // Solve Poisson equation
  std::cout << "Setup function space\n";
  LagrangeElement element(3);
  FESpace Space(mesh, element);

  std::cout << "Assemble matrices and vectors\n";
  BilinearForm Laplace(Space, Space);
  Laplace.AddLaplaceTerm();
  Laplace.AddReactionTerm(c);
  Laplace.Assemble();

  SparseMatrix& Matrix = Laplace.SystemMatrix();
  
  LinearForm F(Space);
  F.AddVolumeForce(f);
  F.Assemble();
  
  Vector& Vec = F.LoadVector();

  std::cout << "Solve equation system\n";
  Vector X(Matrix.Solve(Vec));  

  FEFunction Sol(Space);
  Sol.CreateFunction(X);
  Sol.WriteVtk("solution.vtk");
  
  return 0;
}
