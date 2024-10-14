#include <iostream>
#include <vector>

#include "mesh/Mesh.h"

using namespace chemfem::mesh;

int main()
{
  Mesh mesh("tests/mesh.dat");

  mesh.WriteVtk("mesh_old.vtk", Vector(mesh.NrNodes()));
  
  std::vector<bool> marker(mesh.NrCells(), false);

  marker[0] = true;
  marker[2] = true;

  mesh.Refine(marker);

  std::cout << mesh;
  
  marker[0] = true;
  marker[4] = true;

  mesh.Refine(marker);
  
  mesh.WriteVtk("mesh.vtk", Vector(mesh.NrNodes()));
  return 0;
}
