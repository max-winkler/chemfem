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

      /**
       * Computes the distance to another node
       */
      double Dist(const Node&) const;
      
      friend std::ostream& operator<< (std::ostream&, const Node&);
    private:
      double x, y;
      /* 
       * TODO: It is actually not a good idea to store the index (used also as global DOF)
       * in the Node. The numbering depends on the used FE Space. A mapping between 
       * index and node should be stored there.
       */
      size_t index;
    };
    
  };
};

#endif
