#ifndef _EDGE_H_
#define _EDGE_H_

#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    class Edge
    {
    public:
      Edge(const Cell&, const& Cell);
      
    private:
      const Cell *Neigh0, *Neigh1;
      Node *Node0 = NULL, *Node1 = NULL;
    };
    
  };
};

#endif
