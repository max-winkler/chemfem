#ifndef _EDGE_H_
#define _EDGE_H_

#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    enum EdgeType {BOUNDARY_EDGE, INTERFACE_EDGE};
    
    /**
     * Each instance of this class represents an interface between two cells 
     * of the finite element mesh. Moreover, the endpoints of the edge are stored.
     */
    class Edge
    {           
    public:     
      
      /**
       * Constructor which initializes an edge by the two cells which meet in the edge.
       */
      Edge(Cell&, Cell&);

      /**
       * Constructor which initializes the edge by one cell and the two endpoints.
       * The neighbor of the edge should be set later.
       */
      Edge(Cell&, Node&, Node&);
      
      /**
       * The operator< is used in the mesh class as all edges are stored in an sorted
       * set. An edge is less than another edge if the indices of the nodes fulfill such a relation.
       */
      bool operator<(const Edge&) const;

      /**
       * Updates the neighbor. This method is used when the edge was initialized with 1 Cell only.
       * The neighbor is sometimes found later when iterating over all cells.
       */
      void SetNeighbor(Cell&) const;

      /**
       * Returns a reference to the cells meeting in the edge.
       */
      Cell& GetNeighbor(int) const;

      /**
       * Returns the edge type.
       */
      EdgeType Type() const;
      
      const Node& GetNode(int) const;
    private:
      Node *Node0 = NULL, *Node1 = NULL;
      
      mutable Cell *Neigh0, *Neigh1;
      mutable EdgeType type;      
    };
    
  };
};

#endif
