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
  Mesh mesh("tests/mesh.dat");
  //Mesh mesh = UnitSquareMesh(3);

  const double mu = 0.2;
  
  mesh.WriteVtk("mesh_old.vtk");
  for(int lvl=0; lvl<4; ++lvl)
    mesh.Refine();
  
  for(int lvl=0; lvl<6; ++lvl)
    {
      size_t nr_cells = mesh.NrCells();
      std::cout << "Refinement level " << lvl << ": " << nr_cells << " cells.\n";
      
      std::vector<bool> marker(nr_cells, false);
        
      for(size_t i=0; i<nr_cells; ++i)
        {
	const Cell cell = mesh.GetCell(i);
	const Node center = cell.Barycenter();
	const double h = mesh.MaxDiameter();
	
	double dist = sqrt(pow(center.getX(), 2.) + pow(center.getY(), 2.));

	if(cell.Diameter() > h*pow(dist, 1-mu))
	  marker[i] = true;
        }
      
      mesh.Refine(marker);
    }
  

  if(!mesh.Check())
    {
      std::cerr << "The mesh is broken!\n";
      return -1;
    }
        
  mesh.WriteVtk("mesh.vtk");

  
 
  return 0;
}
  
