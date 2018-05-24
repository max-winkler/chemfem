#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <string>

#include "mesh/Node.h"
#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    /**
     * Mesh class which stores the geometry and the triangulation of the domain
     */
    class Mesh
    {
      friend class UnitSquareMesh;
      
    public:
      /// Returns the number of cells
      size_t NrCells();
      /// Returns the number of nodes
      size_t NrNodes();

      void WriteVtk(const std::string&);
      
    private:
      std::vector<Cell> Cells;
      std::vector<Node> Nodes;
    };
    
  };
};

#endif
