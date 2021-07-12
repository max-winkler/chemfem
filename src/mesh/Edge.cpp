#include "mesh/Edge.h"

namespace chemfem{
  namespace mesh{
    
    Edge::Edge(size_t Node0, size_t Node1, EdgeType type) : Node0(Node0), Node1(Node1), Neigh0(-1), Neigh1(-1), type(type)
    {
    }

    std::pair<size_t, size_t> Edge::GetNeighbors() const
    {
      return std::pair<size_t, size_t>(Neigh0, Neigh1);
    }

    size_t Edge::GetNeighbor(size_t cell) const
    {
      if(Neigh0 == cell)
        return Neigh1;
      else if(Neigh1 == cell)
        return Neigh0;
      else
        std::cerr << "Error in GetNeighbor: The specified cell is not associated with the edge.\n"; 
    }

    void Edge::SetNeighbor(size_t neigh)
    {
      if(Neigh0 == -1)
        Neigh0 = neigh;
      else if(Neigh1 == -1)
        {
	Neigh1 = neigh;
	type = BOUNDARY_EDGE;
        }
      else
        std::cerr << "Error: You try to assign a neighbor to the edge where the neighbors are already set.\n";
    }
    
    EdgeType Edge::Type() const
    {
      return type;
    }

    bool Edge::operator<(const Edge& other) const
    {
      if(Node0 < other.Node0)
	return true;
      if(Node0 == other.Node0 && Node1 < other.Node1)
	return true;
      
      return false;
    }

    bool Edge::operator==(const Edge& other) const
    {
      if(Node0 == other.Node0 && Node1 == other.Node1)
	return true;
      if(Node0 == other.Node1 && Node1 == other.Node0)
	return true;
      
      return false;
    }
    
    std::ostream& operator<< (std::ostream& stream, const Edge& edge)
    {
      stream << "[" << edge.Node0 << "," << edge.Node1 << "]";
					       return stream;
    }
    
  }
}
