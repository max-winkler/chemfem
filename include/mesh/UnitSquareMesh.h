#ifndef _UNIT_SQUARE_MESH_H_
#define _UNIT_SQUARE_MESH_H_

#include <iostream>

#include "mesh/LShapeMesh.h"

namespace chemfem{
  namespace mesh{

    /**
     * Class inherited from Mesh which creates a simple triangulation of the unit square.
     */ 
    class UnitSquareMesh : public Mesh      
    {
    public:
      /**
       * Constructor which creates a mesh for a unit square. The argument is the number of 
       * nodes in each coordinate direction.
       */
      UnitSquareMesh(size_t n);      
    };

  };
};

#endif
