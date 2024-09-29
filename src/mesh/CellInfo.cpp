#include <iomanip>
#include <cmath>

#include "mesh/CellInfo.h"
#include "linalg/Vector.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace mesh{

    CellInfo::CellInfo(const Node& n1, const Node& n2, const Node& n3)
      : LocNode{n1, n2, n3}
    {
    }
    
    double CellInfo::Diam() const
    {
      return std::max(LocNode[1].Dist(LocNode[2]),
                      std::max(LocNode[0].Dist(LocNode[1]),
                               LocNode[0].Dist(LocNode[2]))
                      );
    }

    double CellInfo::Volume() const
    {
      double x[3], y[3];
      GetCellCoords(x, y);
  
      return 0.5*((x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]));
    }

    chemfem::linalg::Vector CellInfo::Barycenter() const
    {
      double x[3], y[3];
      GetCellCoords(x, y);
  
      chemfem::linalg::Vector b(2);
      b[0] = (x[0] + x[1] + x[2])/3.;
      b[1] = (y[0] + y[1] + y[2])/3.;

      return b;
    }

    void CellInfo::GetCellCoords(double* x, double* y) const
    {
      for(int k=0; k<3; ++k)
        {
          x[k] = LocNode[k].getX();
          y[k] = LocNode[k].getY();
        }
    }

    chemfem::linalg::Vector CellInfo::Normal(int k) const
    {
      Vector n(2);

      n[0] = LocNode[(k+1)%3].getY() - LocNode[k].getY();
      n[1] = LocNode[k].getX() - LocNode[(k+1)%3].getX();

      n *= (1./n.Norm());
      return n;
    }

    double CellInfo::EdgeLength(int k) const
    {
      return sqrt(pow(LocNode[(k+1)%3].getX() - LocNode[k].getX(), 2.) +
	        pow(LocNode[(k+1)%3].getY() - LocNode[k].getY(), 2.));
    }
    
    std::ostream& operator<<(std::ostream& os, const CellInfo& info)
    {
      os << "( ";
      for(size_t i=0; i<3; ++i)
        {
          os << info.LocNode[i];
          if(i<2)
            os << " - ";
        }
      return os;
    }
  }
}
