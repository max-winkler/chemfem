#include "mesh/Edge.h"

namespace chemfem{
  namespace mesh{

    Edge::Edge(Node& Node0, Node& Node1) : Node0(&Node0), Node1(&Node1),
					   Neigh0(NULL), Neigh1(NULL),
					   type(INTERFACE_EDGE)
    {
      if(this->Node0->Index() > this->Node1->Index())
	std::swap(this->Node0, this->Node1);
    }    
    
    Edge::Edge(Cell& Neigh0, Cell& Neigh1) : Neigh0(&Neigh0), Neigh1(&Neigh1), type(INTERFACE_EDGE)
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
      : Node0(&Node0), Node1(&Node1), Neigh0(&Neigh0), Neigh1(NULL), type(BOUNDARY_EDGE)
    {
      if(this->Node0->Index() > this->Node1->Index())
	std::swap(this->Node0, this->Node1);
    }

    void Edge::SetNeighbor(Cell& other) const
    {
      if(Neigh1 != NULL)
	std::cerr << "All neighbors of the edge are already set. This seems to be a mistake.\n";

      if(Neigh0 == NULL)
	Neigh0 = &other;
      else
	Neigh1 = &other;
      
      type = INTERFACE_EDGE;
    }

    Cell& Edge::GetNeighbor(int i) const
    {

      if(i!= 0 && i != 1)
	std::cerr << "Requested " << i << "th neighbor of an edge. This makes no sense.\n";

      if(i==0)
	return *Neigh0;
      else 
	return *Neigh1; 
    }

    const Node& Edge::GetNode(int i) const
    {
      return i==0 ? *Node0 : *Node1;
    }

    EdgeType Edge::Type() const
    {
      if(Neigh0 != NULL && Neigh1 != NULL)
	return INTERFACE_EDGE;
      else
	return BOUNDARY_EDGE;
    }

    bool Edge::operator<(const Edge& other) const
    {
      if(Node0->Index() < other.Node0->Index())
	return true;
      if(Node0->Index() == other.Node0->Index() && Node1->Index() < other.Node1->Index())
	return true;
      
      return false;
    }

    
    std::ostream& operator<< (std::ostream& stream, const Edge& edge)
    {
      stream << "[" << *(edge.Node0) << "," << *(edge.Node1) << "]";
					       return stream;
    }
    
  }
}
