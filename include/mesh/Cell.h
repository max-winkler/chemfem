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
     * Class which represents a single cell (triangle/tetrahedron) of a finite element triangulation.
     * This class stored only the indices of the edges and nodes of the cell. On-the-fly, one can create
     * the more detailled Cell class using Mesh::GetCell(size_t).
     */
    class CellInfo
    {
      friend class Edge;
      friend class Mesh;
      friend class chemfem::fem::FESpace;
      friend class chemfem::fem::BilinearForm;
      friend class chemfem::fem::LinearForm;
      friend class chemfem::fem::ErrorNorm;
      
    public:
      /**
       * Constructor creating a new cell by its given vertices.
       */
      CellInfo(size_t, size_t, size_t);

      /// Copy constructor
      CellInfo(const CellInfo&);
      
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
      friend std::ostream& operator<<(std::ostream&, const CellInfo&);
      
    private:
      /// Stores the indices of the 3 vertices of the triangle
      size_t LocNode[3];
      /// Store the indices of the 3 edges of the triangle
      int LocEdge[3];
      /// Stores the index of the cell. Merely used as temporary variable during FESpace::FESpace()
      size_t index;
    };

    /**
     * Triangular cell object storing connectivity information and instances to nodes and edges.
     * Class provides some useful methods used in assembly and error estimation routines.
     */
    class Cell : public CellInfo
    {
    public:
      /// Create a Cell instance
      Cell(const CellInfo&, const Node&, const Node&, const Node&);
      
      /// Returns the coordinates of the element barycenter. Return value is an object of type Node.
      const Node Barycenter() const;

      /// Returns the volume of the cell
      double Volume() const;

      /// Returns the diameter (length of longest edge) of the cell
      double Diameter() const;
      
    private:
      const Node& Node0;
      const Node& Node1;
      const Node& Node2;
    };    
  };
};

#endif
