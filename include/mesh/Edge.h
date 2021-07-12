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
       * Constructor which initializes an egde by two endpoints.
       */
      Edge(Node&, Node&);
       
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
       * Updates the neighbor. This method is used when the edge was initialized with 1 Cell only.
       * The neighbor is sometimes found later when iterating over all cells.
       */
      void SetNeighbor(Cell&) const;

      /**
       * Returns a reference to the cells meeting in the edge.
       */
      Cell& GetNeighbor(int) const;

      /**
       * Returns a reference to the neighbor cell which is not the cell specified in the argument
       */
      Cell& GetNeighbor(const Cell&) const;
      
      /**
       * Returns the edge type.
       */
      EdgeType Type() const;
      
      /**
       * Returns the i-th endpoint of the edge.
       */
      const Node& GetNode(int) const;

      /**
       * The operator< is used in the mesh class as all edges are stored in an sorted
       * set. An edge is less than another edge if the indices of the nodes fulfill such a relation.
       */
      bool operator<(const Edge&) const;

      friend std::ostream& operator<< (std::ostream&, const Edge&);
    private:
      Node *Node0 = NULL, *Node1 = NULL;
      
      mutable Cell *Neigh0, *Neigh1;
      mutable EdgeType type;      
    };
    
  };
};

#endif
