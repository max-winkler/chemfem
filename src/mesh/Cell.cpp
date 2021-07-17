#include <iomanip>
#include <cmath>

#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    CellInfo::CellInfo(size_t n1, size_t n2, size_t n3)
    {
      LocNode[0] = n1;
      LocNode[1] = n2;
      LocNode[2] = n3;
    }

    CellInfo::CellInfo(const CellInfo& other)
      : index(other.index)
    {
      for(int k=0; k<3; ++k)
	{
	  LocNode[k] = other.LocNode[k];
	  LocEdge[k] = other.LocEdge[k];
	}
    }

    std::ostream& operator<<(std::ostream& os, const CellInfo& cell)
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
        
    size_t CellInfo::Index() const
    {
      return index;
    }

    void CellInfo::SetIndex(size_t index)
    {
      this->index = index;
    }

    int CellInfo::EdgeIndex(const size_t edge) const
    {
      size_t edge_ind;
      for(edge_ind=0; edge_ind<3; ++edge_ind)
	if(LocEdge[edge_ind] == edge) break;

      if(edge_ind == 3) edge_ind = -1;

      return edge_ind;
    }

    // Methods for class Cell
    
    Cell::Cell(const CellInfo& info, const Node& n0, const Node& n1, const Node& n2)
      : CellInfo(info), Node0(n0), Node1(n1), Node2(n2)
    {      
    }

    const Node Cell::Barycenter() const
    {
      double x = (Node0.getX() + Node1.getX() + Node2.getX())/3.;
      double y = (Node0.getY() + Node1.getY() + Node2.getY())/3.;

      return Node(0, x, y);
    }

    double Cell::Volume() const
    {
      double x0 = Node1.getX() - Node0.getX();
      double x1 = Node2.getX() - Node0.getX();
      double y0 = Node1.getY() - Node0.getY();
      double y1 = Node2.getY() - Node0.getY();

      return 0.5*(x0*y1-y0*x1);
    }

    double Cell::Diameter() const
    {
      double x0 = Node1.getX() - Node0.getX();
      double y0 = Node1.getY() - Node0.getY();
      double x1 = Node2.getX() - Node0.getX();
      double y1 = Node2.getY() - Node0.getY();
      double x2 = Node2.getX() - Node1.getX();
      double y2 = Node2.getY() - Node1.getY();
      
      double l0 = sqrt(x0*x0+y0*y0);
      double l1 = sqrt(x1*x1+y1*y1);
      double l2 = sqrt(x2*x2+y2*y2);

      return std::max(l0, std::max(l1, l2));
    }

  };
};
