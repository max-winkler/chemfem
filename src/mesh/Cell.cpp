#include "mesh/Cell.h"

namespace chemfem{
  namespace mesh{

    Cell::Cell(Node& n1, Node& n2, Node& n3)
    {
      LocNode[0] = &n1;
      LocNode[1] = &n2;
      LocNode[2] = &n3;
    }

    double Cell::Determinant() const
    {
      double x[3], y[3];
      for(int k=0; k<3; ++k)
	{
	  x[k] = LocNode[k]->getX();
	  y[k] = LocNode[k]->getY();
	}
      
      return (x[1] - x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]);
    }
  };
};
