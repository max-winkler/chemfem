#ifndef _NODE_H_
#define _NODE_H_

#include <iostream>

namespace chemfem{
  namespace fem{
    class FESpace;
  }
}

namespace chemfem{
  namespace mesh{

    /**
     * Class represents a single node of a finite element triangulation.
     */
    class Node
    {
      friend class Mesh;
      friend class chemfem::fem::FESpace;

    public:
      /**
       * Constructor which initializes the node by its coordinates.
       */
      Node(size_t, double, double);

      /**
       * Returns the x-coordinate
       */
      double getX() const;

      /**
       * Returns the y-coordinate
       */
      double getY() const;

      /**
       * Returns the index of the node.
       */
      size_t Index() const;
    private:
      double x, y;
      size_t index;
    };
    
  };
};

#endif
