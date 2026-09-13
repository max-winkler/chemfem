#include <iostream>
#include <vector>

#include "mesh/Mesh.h"

using namespace chemfem::mesh;

// Reads a mesh in the gmsh format and refines it locally. The file is written in the
// format of gmsh 2.2, the meshes of a current gmsh use version 4.1.

int main()
{
  Mesh mesh("meshes/simple.msh");

  std::cout << mesh.NrNodes() << " nodes, " << mesh.NrCells() << " cells, "
            << mesh.NrEdges() << " edges" << std::endl;

  bool ok = true;

  if(mesh.NrNodes() != 5 || mesh.NrCells() != 3)
    {
      std::cerr << "ERROR: the mesh file was not read correctly.\n";
      ok = false;
    }

  if(!mesh.Check())
    {
      std::cerr << "ERROR: the mesh read from the file is broken.\n";
      ok = false;
    }

  for(int level=0; level<3 && ok; ++level)
    {
      std::vector<bool> marker(mesh.NrCells(), false);
      marker[0] = true;
      marker[mesh.NrCells()/2] = true;

      mesh.Refine(marker);

      std::cout << "refinement " << level+1 << ": " << mesh.NrCells() << " cells"
                << std::endl;

      if(!mesh.Check())
        {
          std::cerr << "ERROR: the mesh is broken after refinement " << level+1 << ".\n";
          ok = false;
        }
    }

  if(!ok)
    return 1;

  std::cout << "\nMeshTest was successful.\n";
  return 0;
}
