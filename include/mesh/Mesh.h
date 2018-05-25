#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <string>

#include "mesh/Node.h"
#include "mesh/Cell.h"

// Forward declarations of friend classes
namespace chemfem::fem{
  class FESpace;
  class BilinearForm;
}

namespace chemfem{
  namespace mesh{

    /**
     * Mesh class which stores the geometry and the triangulation of the domain
     */
    class Mesh
    {
      friend class UnitSquareMesh;
      friend class chemfem::fem::FESpace;
      friend class chemfem::fem::BilinearForm;
    public:
      /// Returns the number of cells
      size_t NrCells() const;
      /// Returns the number of nodes
      size_t NrNodes() const;

      void WriteVtk(const std::string&);
      
    private:
      std::vector<Cell> Cells;
      std::vector<Node> Nodes;
    };
    
  };
};

#endif
