#ifndef _LSHAPE_MESH_H_
#define _LSHAPE_MESH_H_

#include <iostream>

#include "mesh/Mesh.h"

namespace chemfem{
  namespace mesh{

    /**
     * Class inherited from Mesh which creates a simple triangulation of an L-shaped domain.
     */
    class LShapeMesh : public Mesh
    {
    public:
      /**
       * Constructor creating a mesh for an L-shaped domain. The argument is the number
       * of nodes in each coordinate direction.
       */
      LShapeMesh(size_t n);
    };

  };
};

#endif
