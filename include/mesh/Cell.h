#ifndef _CELL_H_
#define _CELL_H_

#include <iostream>

#include "Node.h"

// Forward declarations for friend classes
namespace chemfem::fem{
  class FESpace;
}

namespace chemfem{
  namespace mesh{

    /**
     * Class which represents a single cell (triangle/tetrahedron) of a finite element triangulation.
     */
    class Cell
    {
      friend class Mesh;
      friend class chemfem::fem::FESpace;
      
    public:
      /**
       * Constructor creating a new cell by its given vertices.
       */
      Cell(Node&, Node&, Node&);

      /**
       * Returns the volume of the parallelogram which is spanned by the vertices of the cell.
       */
      double Determinant() const;
     
    private:
      /// Stores the 3 vertices of the triangle
      Node *LocNode[3];
    };
    
  };
};

#endif
