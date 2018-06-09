#include "mesh/Edge.h"

namespace chemfem{
  namespace mesh{

    Edge::Edge(const Cell& Neigh0, const Cell& Neigh0) : Neigh0(&Neigh0), Neigh1(&Neigh1)
    {
      for(Node* Node0 = Neigh0.LocNode; Node0 != Neigh0.LocNode+3; ++Node0)
	for(Node* Node1 = Neigh1.LocNode; Node1 != Neigh1.LocNode+3; ++Node1)
	  {
	    if(Node0 == Node1)
	      {
		if(this->Node0 == NULL)
		  this->Node0 = Node0;
		else
		  this->Node1 = Node0;
	      }	    
	  }
      if(this->Node1 == NULL)
	std::cerr << "Broken mesh. Can not find joint vertices of the two Cells.\n";
    }    
    
  }
}
