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
  // Mesh refinement parameters
  double mu     = 0.1;     // refinement strength
  double d      = 0.5;     // refinement radius
  int levels    = 10;      // maximum refinement levels
    
  // Build up and refine mesh
  LShapeMesh mesh(4);
  
  mesh.WriteVtk("mesh_old.vtk");

  for(int i=0; i<levels; ++i)
    {
      std::cout << "Refinement level " << i+1 << " of " << levels << std::endl;
      auto Cells = mesh.GetCellList();
      
      std::vector<bool> marker(mesh.NrCells(), false);

      // Compute global mesh size
      // \todo Provide method in Mesh class to do this
      double h = 0;
      for(size_t c=0; c<mesh.NrCells(); ++c)
        {
          CellInfo info = mesh.GetCellInfo(c);

          double hT = info.Diam();

          if(hT > h)
            h = hT;
        }

      // Mark elements to be refined according to the mesh refinement criterion
      for(size_t c=0; c<mesh.NrCells(); ++c)
        {
          CellInfo info = mesh.GetCellInfo(c);

          double hT = info.Diam();
          double rT = (info.Barycenter() - Coordinate{0., 0.}).Norm();
          
          if(hT > pow(rT/d, 1.-mu) * h)
            {
              std::cout << "  cell " << c << " has been marked for refinement.\n";
              marker[c] = true;
            }
        }
      
      mesh.Refine(marker);
      if(!mesh.Check())
        {
          std::cerr << "The mesh is broken!\n";
          return -1;
        }
    }
  mesh.WriteVtk("mesh.vtk");
  
  return 0;
}
  
