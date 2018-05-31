#include <iostream>
#include <string>

#include "quadrature/QuadFormula.h"

namespace chemfem{
  namespace quadrature{

    QuadratureFormula::QuadratureFormula(QUAD_FORMULA formula)
    {
      switch(formula)
	{
	case MIDPOINT:
	  points = 1;
	  weights = new double[1];
	  xi = new double[1];
	  eta = new double[1];

	  xi[0] = 1./3;
	  eta[0] = 1./3;
	  weights[0] = 1.;
	  
	  break;
	default:
	  std::cerr << "Quadrature formula not implemented yet\n";
	}
    }

    QuadratureFormula::QuadratureFormula(int degree)
    {
      std::cerr << "Initialization of quadrature formula by degree not implemented yet\n";
    }

    size_t QuadratureFormula::NrQuadPoints()
    {
      return points;
    }
    
    void QuadratureFormula::FormulaData(double*& weights_in, double*& xi_in, double*& eta_in)
    {
      weights_in = weights;
      xi_in = xi;
      eta_in = eta;
    }

  }
}
