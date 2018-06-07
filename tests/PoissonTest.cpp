#include <iostream>

#include "fem/BilinearForm.h"
#include "fem/LagrangeElement.h"
#include "mesh/UnitSquareMesh.h"
#include "linalg/SparseMatrix.h"

using namespace chemfem::fem;
using namespace chemfem::linalg;
using namespace chemfem::mesh;

int main()
{
  UnitSquareMesh mesh(3);

  std::cout << "Nr of nodes : " << mesh.NrNodes() << std::endl;
  std::cout << "Nr of cells : " << mesh.NrCells() << std::endl;
  
  LagrangeElement element(1);
  FESpace Space(mesh, element);

  BilinearForm Laplace(Space, Space);
  Laplace.AddLaplaceTerm();
  Laplace.Assemble();

  SparseMatrix& Matrix = Laplace.SystemMatrix();
  
  std::cout << "The finite element system matrix is:\n" << Matrix << std::endl;

  return 0;
}
