#include <cstring>

#include "mesh/Edge.h"

namespace chemfem{
  namespace mesh{

    Edge::Edge() : Node0(0), Node1(1), Neigh0(-1), Neigh1(-1), type(BOUNDARY_EDGE) {}
    
    Edge::Edge(size_t Node0, size_t Node1, EdgeType type) : Node0(Node0), Node1(Node1), Neigh0(-1), Neigh1(-1), type(type)
    {
    }

    std::pair<size_t, size_t> Edge::GetNeighbors() const
    {
      return std::pair<size_t, size_t>(Neigh0, Neigh1);
    }

    size_t Edge::GetNeighbor(size_t cell) const
    {
      if(Neigh0 == (long)cell)
        return Neigh1;
      else if(Neigh1 == (long)cell)
        return Neigh0;

      std::cerr << "Error in GetNeighbor: The specified cell with index " << cell
	      << " is not associated with the edge. The edge connects " << Neigh0 << " and "
	      << Neigh1 << ".\n";
      return -1;
    }

    size_t Edge::GetNeighborByLocalIndex(int cell) const
    {
      if(cell == 0)
        return Neigh0;
      else if(cell == 1)
        return Neigh1;

      std::cerr << "Error in GetNeighborByIndex: The specified cell with local index " << cell
		<< " is not associated with the edge.\n";
      return -1;
    }

    void Edge::SetNeighbor(size_t neigh)
    {
      const long n = (long)neigh;

      if(Neigh0 == n || Neigh1 == n)
        return;

      if(Neigh0 == -1)
        Neigh0 = n;
      else if(Neigh1 == -1)
        {
	Neigh1 = n;
	type = INTERFACE_EDGE;
        }
      else
        {
	std::cerr << "Error: You try to assign a neighbor to the edge where the neighbors are already set.\n";
	std::cerr << "  Node0  = " << Node0  << ", Node1 = "  << Node1  << std::endl;
	std::cerr << "  Neigh0 = " << Neigh0 << ", Neigh1 = " << Neigh1 << ", neigh = " << neigh << std::endl;
        }
    }

    void Edge::UnsetNeighbor(size_t neigh)
    {
      if(Neigh0 == (long)neigh)
        Neigh0 = -1;
      else if(Neigh1 == (long)neigh)
        Neigh1 = -1;
      else
        std::cerr << "Error: The cell with index " << neigh << " is not adjacent to the edge.\n";
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
