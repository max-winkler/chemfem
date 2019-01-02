#include <iostream>

#include "mesh/UnitSquareMesh.h"

using namespace chemfem::mesh;

int main()
{
  Mesh mesh = UnitSquareMesh(10);
  mesh.WriteVtk("mesh_old.vtk");
  mesh = mesh.RefineUniform();
  
  mesh.WriteVtk("mesh.vtk");
  return 0;
}
