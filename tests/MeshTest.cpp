#include <iostream>

#include "mesh/UnitSquareMesh.h"

using namespace chemfem::mesh;

int main()
{
  Mesh mesh(UnitSquareMesh(10));
  mesh.WriteVtk("mesh.vtk", Vector(mesh.NrNodes()));
  return 0;
}
