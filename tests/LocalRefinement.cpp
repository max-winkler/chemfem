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

int main()
{
  // Build up and refine mesh
  Mesh mesh = UnitSquareMesh(4);

  mesh.WriteVtk("mesh_old.vtk");
  std::cout << mesh << std::endl;
 
  std::vector<bool> marker(mesh.NrCells(), false);

  marker[2] = true;
  marker[3] = true;
  marker[4] = true;
  marker[8] = true;
  marker[9] = true;
 
  mesh = mesh.Refine(marker);
  if(!mesh.Check())
    {
      std::cerr << "The mesh is broken!\n";
      return -1;
    }
        
  mesh.WriteVtk("mesh.vtk");
  
  return 0;
}
  
