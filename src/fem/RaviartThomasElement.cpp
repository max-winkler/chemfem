#include "fem/RaviartThomasElement.h"

#include <iostream>

using chemfem::linalg::Coordinate;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    RaviartThomasElement::RaviartThomasElement() : VectorElement(FEType::RaviartThomas, 0)
    {
      nr_dof = 3;
      dofs_per_edge = 1;
    }

    Vector2D RaviartThomasElement::Value(int k, double xi, double eta) const
    {
      switch(k)
	{
	case 0: return Vector2D(xi, eta - 1.);
	case 1: return Vector2D(xi, eta);
	case 2: return Vector2D(xi - 1., eta);
	}

      std::cerr << "Requested basis function " << k
		<< ", but the element has only " << nr_dof << " degrees of freedom.\n";

      return Vector2D(0., 0.);
    }

    Matrix2D RaviartThomasElement::Gradient(int k, double, double) const
    {
      if(k < 0 || k >= nr_dof)
	{
	  std::cerr << "Requested basis function " << k
		    << ", but the element has only " << nr_dof << " degrees of freedom.\n";
	  return Matrix2D(0., 0., 0., 0.);
	}

      return Matrix2D(1., 0., 0., 1.);
    }

    Coordinate RaviartThomasElement::NodalPoint(int k) const
    {
      double xi, eta;
      EdgeToRefCoords(k, 0.5, xi, eta);

      return Coordinate{xi, eta};
    }

  }
}
