#include <iomanip>
#include <cmath>

#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    Cell::Cell(size_t n1, size_t n2, size_t n3)
    {
      LocNode[0] = n1;
      LocNode[1] = n2;
      LocNode[2] = n3;
    }

    Cell::Cell(const Cell& other)
      : index(other.index)
    {
      for(int k=0; k<3; ++k)
	{
	  LocNode[k] = other.LocNode[k];
	  LocEdge[k] = other.LocEdge[k];
	}
    }

    std::ostream& operator<<(std::ostream& os, const Cell& cell)
    {
      os << "( ";
      for(size_t i=0; i<3; ++i)
        {
	os << std::setw(3) << cell.LocNode[i];
	if(i<2)
	  os << " , ";
        }
      os << " ) - edges: ( ";
      for(size_t i=0; i<3; ++i)
	{
	  os << std::setw(3) << cell.LocEdge[i];
	  if(i<2)
	    os << " , ";
	}
      os << " )";
      return os;
    }
        
    size_t Cell::Index() const
    {
      return index;
    }

    void Cell::SetIndex(size_t index)
    {
      this->index = index;
    }

    int Cell::EdgeIndex(const size_t edge) const
    {
      size_t edge_ind;
      for(edge_ind=0; edge_ind<3; ++edge_ind)
	if(LocEdge[edge_ind] == edge) break;

      if(edge_ind == 3) edge_ind = -1;

      return edge_ind;
    }
        
  };
};
