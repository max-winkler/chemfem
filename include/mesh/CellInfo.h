#ifndef _CELL_INFO_H_
#define _CELL_INFO_H_

#include <iostream>

#include "mesh/Node.h"

#include "linalg/Coordinate.h"
#include "linalg/Vector2D.h"

namespace chemfem{
  namespace mesh{

    /**
     * Structure storing information on a mesh cell and providing auxiliary functions for 
     * the computation of geometric quantities like cell diameter, volume and barycenter.
     * A CellInfo object is usually generated on-the-fly by Mesh::GetCellInfo(size_t).
     */
    class CellInfo
    {
      friend class Mesh;
      
    public:

      /**
       * Create CellInfo object by specifying the vertices of the cell
       */
      CellInfo(const Node&, const Node&, const Node&);
      
      /**
       * Computes the diameter of the cell, meaning the maximum distance of two nodes.
       */ 
      double Diam() const;

      /**
       * Computes the surface area of the cell.
       */
      double Volume() const;

      /**
       * Returns the element's barycenter
       */
      chemfem::linalg::Coordinate Barycenter() const;

      /**
       * Returns the normal vector at the k-th edge
       */
      chemfem::linalg::Vector2D Normal(int) const;

      /**
       * Returns length of the k-th edge
       */
      double EdgeLength(int) const;
      
      /**
       * Console output of the cell info
       */
      friend std::ostream& operator<<(std::ostream&, const CellInfo&);
      
    private:
      Node LocNode[3];

      /**
       * Helper function to store node coordinates in vectors
       */
      void GetCellCoords(double*, double*) const;       
    };
    
  }
}

#endif
