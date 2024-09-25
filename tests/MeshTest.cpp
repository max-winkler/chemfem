#include <iostream>

#include "mesh/UnitSquareMesh.h"

using namespace chemfem::mesh;

int main()
{
  LShapeMesh mesh(5);

  mesh.Refine();
  mesh.Refine();
  
  mesh.WriteVtk("mesh.vtk", Vector(mesh.NrNodes()));
  return 0;
}
