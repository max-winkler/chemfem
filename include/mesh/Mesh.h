#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <set>
#include <string>

#include "mesh/Node.h"
#include "mesh/Cell.h"
#include "mesh/Edge.h"

#include "linalg/Vector.h"

using chemfem::linalg::Vector;

// Forward declarations of friend classes
namespace chemfem::fem{
  class FESpace;
  class BilinearForm;
  class LinearForm;
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
      friend class chemfem::fem::LinearForm;
    public:

      /**
       * Destructor which frees all allocated memory by cells, nodes and edges.
       */
      ~Mesh();
      
      /// Returns the number of cells
      size_t NrCells() const;
      /// Returns the number of nodes
      size_t NrNodes() const;
      /// Write Mesh to a VTK file
      void WriteVtk(const std::string&, const Vector& x) const;

      /**
       * Refine the mesh uniformly (all triangles are decomposed into 4 of equivalent size)
       */
      void RefineUniform();	

      /**
       * Refine the mesh according to the refinement description of the cells.
       */
      void Refine();
      
      /**
       * Returns a reference to the cell list of the mesh.
       */
      const std::vector<Cell>& GetCellList() const;
      
    private:
      std::vector<Cell> Cells;
      std::vector<Node> Nodes;
      std::set<Edge> Edges;
      std::vector<const Edge*> BdEdges;
      
      void CreateEdgeList();
      
    };
    
  };
};

#endif
