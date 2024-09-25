#ifndef _MESH_H_
#define _MESH_H_

#include <iostream>
#include <vector>
#include <set>
#include <string>

#include "mesh/Node.h"
#include "mesh/Cell.h"
#include "mesh/Edge.h"
#include "mesh/CellInfo.h"

#include "linalg/Vector.h"

using chemfem::linalg::Vector;

// Forward declarations of friend classes
namespace chemfem::fem{
  class FESpace;
  class BilinearForm;
  class LinearForm;
  class ErrorNorm;
}

typedef std::pair<std::set<chemfem::mesh::Edge>::iterator, bool> set_insert_res;

namespace chemfem{
  namespace mesh{

    /**
     * Mesh class which stores the geometry and the triangulation of the domain
     */
    class Mesh
    {
      // Friend declarations
      friend class chemfem::fem::FESpace;
      friend class chemfem::fem::BilinearForm;
      friend class chemfem::fem::LinearForm;
      friend class chemfem::fem::ErrorNorm;
    public:

      /// Default constructor
      Mesh();
      
      /// Copy constructor
      Mesh(const Mesh&);

      /// Load mesh from file
      Mesh(const std::string&);
      
      /**
       * Destructor which frees all allocated memory by cells, nodes and edges.
       */
      ~Mesh();

      /// Copy assignment
      Mesh& operator=(const Mesh&);
      
      /// Returns the number of cells
      size_t NrCells() const;
      /// Returns the number of nodes
      size_t NrNodes() const;
      /// Returns the number of edges
      size_t NrEdges() const;
      
      /// Write Mesh to a VTK file
      void WriteVtk(const std::string&) const;      
      /// Write Mesh to a VTK file
      void WriteVtk(const std::string&, const Vector& x) const;

      /**
       * Refine the mesh uniformly (all triangles are decomposed into 4 of equivalent size)
       */
      void RefineUniform();	

      /**
       * Refine all elements
       */
      void Refine();
      
      /**
       * Refine the mesh according to the refinement description of the cells.
       */
      void Refine(const std::vector<bool>&);

      /**
       * Returns a reference to the cell list of the mesh.
       */
      const std::vector<Cell>& GetCellList() const;

      /**
       * Check if the mesh data structure is broken
       */
      bool Check();

      /**
       * Returns the volume of the parallelogram which is spanned by the vertices of the cell.
       */
      double Determinant(size_t) const;

      /**
       * Returns the Jacobian of the reference transformation.
       */
      chemfem::linalg::DenseMatrix Jacobian(size_t) const;

      /**
       * Create CellInfo object providing geometric computations for a mesh cell.
       */
      CellInfo GetCellInfo(size_t) const;
            
      /**
       * Console output of the mesh information
       */
      friend std::ostream& operator<<(std::ostream&, const Mesh&);

    private:

      /// Copies a mesh and updates the pointers to nodes and edges
      void copy(const Mesh&);

    protected:
      void CreateEdgeList();
      
      std::vector<Node> Nodes;
      std::vector<Cell> Cells;
      std::vector<Edge> Edges;

    };
    
  };
};

#endif
