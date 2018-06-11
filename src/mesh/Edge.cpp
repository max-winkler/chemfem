#include "mesh/Edge.h"

namespace chemfem{
  namespace mesh{

    Edge::Edge(Cell& Neigh0, Cell& Neigh1) : Neigh0(&Neigh0), Neigh1(&Neigh1)
    {
      for(Node* Node0 = Neigh0.LocNode[0]; Node0 != Neigh0.LocNode[3]; ++Node0)
	for(Node* Node1 = Neigh1.LocNode[0]; Node1 != Neigh1.LocNode[3]; ++Node1)
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

      if(Node0->Index() > Node1->Index())
	std::swap(Node0, Node1);
    }

    Edge::Edge(Cell& Neigh0, Node& Node0, Node& Node1)
      : Neigh0(&Neigh0), Neigh1(NULL), Node0(&Node0), Node1(&Node1)    
    {
      if(this->Node0->Index() > this->Node1->Index())
	std::swap(this->Node0, this->Node1);
    }
    
    bool Edge::operator==(const Edge& other) const
    {
      return (Node0 == other.Node0 && Node1 == other.Node1)
	|| (Node0 == other.Node1 && Node1 == other.Node0);
    }

    bool Edge::operator<(const Edge& other) const
    {
      if(Node0->Index() < other.Node0->Index())
	return true;
      if(Node0->Index() == other.Node0->Index() && Node1->Index() < other.Node1->Index())
	return true;
      return false;
    }

    void Edge::SetNeighbor(Cell& other) const
    {
      if(Node1 == NULL)
	std::cerr << "All neighbors of the edge are already set. This seems to be a mistake.\n";

      Neigh1 = &other;
    }

    Cell& Edge::GetNeighbor(int i) const
    {
      if(i==0) return *Neigh0;
      else if(i==1) return *Neigh1;
      else
	std::cerr << "Requested " << i << "th neighbor of an edge. This makes no sense.\n";
    }
  }
}
