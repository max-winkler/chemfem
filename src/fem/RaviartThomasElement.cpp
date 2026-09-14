#include "fem/RaviartThomasElement.h"

#include <iostream>

using chemfem::linalg::Coordinate;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    RaviartThomasElement::RaviartThomasElement() : Element(FEType::RaviartThomas, 0)
    {
      nr_dof = 3;
      dofs_per_edge = 1;
      mapping = ContravariantPiola;
    }

    VectorRefValues RaviartThomasElement::VectorReference(int k, double xi, double eta) const
    {
      switch(k)
	{
	case 0: return VectorRefValues{Vector2D(xi, eta - 1.), 2.};
	case 1: return VectorRefValues{Vector2D(xi, eta), 2.};
	case 2: return VectorRefValues{Vector2D(xi - 1., eta), 2.};
	}

      std::cerr << "Requested basis function " << k
		<< ", but the element has only " << nr_dof << " degrees of freedom.\n";

      return VectorRefValues{Vector2D(0., 0.), 0.};
    }

    double RaviartThomasElement::Value(int, double, double) const
    {
      std::cerr << "The Raviart-Thomas basis functions are vector valued, use "
		<< "VectorReference.\n";
      return 0.;
    }

    Vector2D RaviartThomasElement::Gradient(int, double, double) const
    {
      std::cerr << "The gradient of the Raviart-Thomas basis functions is not implemented.\n";
      return Vector2D(0., 0.);
    }

    Matrix2D RaviartThomasElement::Hessian(int, double, double) const
    {
      std::cerr << "The Hessian of the Raviart-Thomas basis functions is not implemented.\n";
      return Matrix2D(0., 0., 0., 0.);
    }

    Coordinate RaviartThomasElement::NodalPoint(int k) const
    {
      double xi, eta;
      EdgeToRefCoords(k, 0.5, xi, eta);

      return Coordinate{xi, eta};
    }

  }
}
