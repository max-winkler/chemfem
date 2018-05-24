#ifndef _CELL_H_
#define _CELL_H_

#include <iostream>

namespace chemfem{
  namespace mesh{

    /**
     * Class which represents a single cell (triangle/tetrahedron) of a finite element triangulation.
     */
    class Cell
    {
      friend class Mesh;
      
    public:
      /**
       * Constructor creating a new cell by its given vertices.
       */
      Cell(size_t, size_t, size_t);
    private:
      /// Stores the 3 vertices of the triangle
      size_t LocIndex[3];
    };
    
  };
};

#endif
