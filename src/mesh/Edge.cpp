#include <cstring>

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
        std::cerr << "Error in GetNeighbor: The specified cell with index " << cell
	        << " is not associated with the edge. The edge connects " << Neigh0 << " and "
	        << Neigh1 << ".\n"; 
    }

    size_t Edge::GetNeighborByLocalIndex(int cell) const
    {
      if(cell == 0)
        return Neigh0;
      else if(cell == 1)
        return Neigh1;
      else
        std::cerr << "Error in GetNeighborByIndex: The specified cell with local index " << cell
		  << " is not associated with the edge.\n"; 
    }

    void Edge::SetNeighbor(size_t neigh)
    {
      if(Neigh0 == neigh)
        return;
      
      if(Neigh0 == -1)
        Neigh0 = neigh;
      else if(Neigh1 == -1 && Neigh1 != neigh)
        {
	Neigh1 = neigh;
	type = INTERFACE_EDGE;
        }
      else if(Neigh1 == neigh || Neigh0 == neigh) {}
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
    
    std::ostream& operator<<(std::ostream& os, const Edge& edge)
    {
      std::string type_str = edge.type==BOUNDARY_EDGE ? "boundary edge" : "interface edge"; 
      os << " nodes: ( " << edge.Node0 << " , " << edge.Node1 << " ),  cells: ( "
	 << edge.Neigh0 << " , " << edge.Neigh1 << " ) - " << type_str;
      return os;
    }    
  }
}
