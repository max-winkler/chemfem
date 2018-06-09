#ifndef _CELL_H_
#define _CELL_H_

#include <iostream>

#include "linalg/DenseMatrix.h"
#include "Node.h"

// Forward declarations for friend classes
namespace chemfem::fem{
  class FESpace;
  class BilinearForm;
  class LinearForm;
}

namespace chemfem{
  namespace mesh{

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
      
    private:
      /// Stores the 3 vertices of the triangle
      Node *LocNode[3];
    };
    
  };
};

#endif
