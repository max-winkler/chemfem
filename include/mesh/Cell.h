#ifndef _CELL_H_
#define _CELL_H_

#include <iostream>

#include "linalg/DenseMatrix.h"

#include "mesh/Node.h"
#include "mesh/Edge.h"

// Forward declarations for friend classes
namespace chemfem::fem{
  class FESpace;
  class BilinearForm;
  class LinearForm;
}

namespace chemfem{
  namespace mesh{

    // Forward declaration to dependent classes
    class Edge;
    
    /**
     * Class which represents a single cell (triangle/tetrahedron) of a finite element triangulation.
     */
    class Cell
    {
      friend class Edge;
      friend class Mesh;
      friend class chemfem::fem::FESpace;
      friend class chemfem::fem::BilinearForm;
      friend class chemfem::fem::LinearForm;
      
    public:
      /**
       * Constructor creating a new cell by its given vertices.
       */
      Cell(Node&, Node&, Node&);

      /**
       * Returns the volume of the parallelogram which is spanned by the vertices of the cell.
       */
      double Determinant() const;

      /**
       * Returns the Jacobian of the reference transformation.
       */
      chemfem::linalg::DenseMatrix Jacobian() const;

      /**
       * Returns the index Variable
       */
      size_t Index() const;

      /**
       * Set the index variable.
       */
      void SetIndex(size_t);

      /**
       * Returns the local index of the edge in the cell. The return value is 
       * -1 if the edge is not an edge of the cell.
       */
      int EdgeIndex(const Edge&) const;
      
    private:
      /// Stores the 3 vertices of the triangle
      Node *LocNode[3];
      /// Store pointers to the 3 edges of the triangle
      const Edge *LocEdge[3];
      
      /// Stores the index of the cell. Merely used as temporary variable during FESpace::FESpace()
      size_t index;
    };
    
  };
};

#endif
