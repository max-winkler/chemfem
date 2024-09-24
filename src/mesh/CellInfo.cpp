#include <iomanip>

#include "mesh/CellInfo.h"

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
