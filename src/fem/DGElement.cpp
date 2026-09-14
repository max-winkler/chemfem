#include "fem/DGElement.h"

#include <iostream>

using chemfem::linalg::Coordinate;
using chemfem::linalg::Matrix2D;
using chemfem::linalg::Vector2D;

namespace chemfem{
  namespace fem{

    DGElement::DGElement(int degree) : Element(FEType::DG, degree)
    {
      switch(degree)
	{
	case 0: nr_dof = 1; break;
	default:
	  nr_dof = 0;
	  std::cerr << "DG elements of degree " << degree << " are not implemented yet\n";
	}

      dofs_interior = nr_dof;
    }

    double DGElement::Value(int, double, double) const
    {
      return 1.;
    }

    Vector2D DGElement::Gradient(int, double, double) const
    {
      return Vector2D(0., 0.);
    }

    Matrix2D DGElement::Hessian(int, double, double) const
    {
      return Matrix2D(0., 0., 0., 0.);
    }

    Coordinate DGElement::NodalPoint(int) const
    {
      return Coordinate{1./3, 1./3};
    }

  }
}
