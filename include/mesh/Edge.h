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
      friend class Mesh;
      friend class chemfem::fem::FESpace;
    public:     

      /**
       * Constructor which initializes an edge by the two cells which meet in the edge.
       */
      Edge(size_t, size_t, EdgeType);
      
      /**
       * Returns a reference to the cells meeting in the edge.
       */
      std::pair<size_t, size_t> GetNeighbors() const;

      /**
       * Returns a reference to the neighbor cell which is not the cell specified in the argument
       */
      size_t GetNeighbor(size_t) const;

      /**
       * Assigns a neighbor to the edge
       */
      void SetNeighbor(size_t);
      
      /**
       * Returns the edge type.
       */
      EdgeType Type() const;
      
      /**
       * The operator< is used in the mesh class as all edges are stored in an sorted
       * set. An edge is less than another edge if the indices of the nodes fulfill such a relation.
       */
      bool operator<(const Edge&) const;

      /**
       * Two edges are equal when the indices of the two nodes are equal
       */
      bool operator==(const Edge&) const;
      
      friend std::ostream& operator<< (std::ostream&, const Edge&);
    private:
      size_t Node0, Node1;
      mutable long Neigh0, Neigh1;
      mutable EdgeType type;      
    };
    
  };
};

#endif
