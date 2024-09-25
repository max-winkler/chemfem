#ifndef _CELL_H_
#define _CELL_H_

#include <iostream>

#include "linalg/DenseMatrix.h"

#include "mesh/Node.h"
#include "mesh/Edge.h"
#include "mesh/RefData.h"

// Forward declarations for friend classes
namespace chemfem::fem{
  class FESpace;
  class BilinearForm;
  class LinearForm;
  class ErrorNorm;
}

namespace chemfem{
  namespace mesh{

    // Forward declaration to dependent classes
    class Edge;
    
    /**
     * Triangular cell object storing connectivity information and instances to nodes and edges.
     * Class provides some useful methods used in assembly and error estimation routines.
     */
    class Cell
    {      
      friend class Edge;
      friend class Mesh;
      friend class CellInfo;
      friend class chemfem::fem::FESpace;
      friend class chemfem::fem::BilinearForm;
      friend class chemfem::fem::LinearForm;
      friend class chemfem::fem::ErrorNorm;
      
    public:
      /**
       * Constructor creating a new cell by its given vertices.
       */
      Cell(size_t, size_t, size_t);

      /// Copy constructor
      Cell(const Cell&);
      
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
      int EdgeIndex(const size_t) const;

      /**
       * Console output of the cell
       */
      friend std::ostream& operator<<(std::ostream&, const Cell&);
      
    private:
      /// Stores the indices of the 3 vertices of the triangle
      size_t LocNode[3];
      /// Store the indices of the 3 edges of the triangle
      int LocEdge[3];
      /// Stores the index of the cell. Merely used as temporary variable during FESpace::FESpace()
      size_t index;
      /// A flag to figure our if the cell was obtained from the bisection of another cell
      size_t orig;
    };
    
  };
};

#endif
