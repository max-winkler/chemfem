#include <iostream>
#include <string>

#include "quadrature/QuadFormula.h"

using chemfem::linalg::Vector;

namespace chemfem{
  namespace quadrature{

    QuadratureFormula::QuadratureFormula(QUAD_FORMULA formula)
    {
      switch(formula)
	{
	case MIDPOINT:
	  Points = 1;
	  Weights = Vector(1);
	  Xi = Vector(1);
	  Eta = Vector(1);

	  Xi[0] = 1./3;
	  Eta[0] = 1./3;
	  Weights[0] = 1.;
	  
	  break;
	case VERTEX:
	  Points = 3;

	  Weights = Vector(3);
	  Xi = Vector(3);
	  Eta = Vector(3);

	  for(int i=0; i<Points; ++i)
	    Weights[i] = 1./3;
	  
	  Xi[0] = 0.; Eta[0] = 0.;
	  Xi[1] = 1.; Eta[1] = 0.;
	  Xi[2] = 0.; Eta[2] = 1.;  	    

	  break;

	case GAUSS_EDGE:
	  Points = 3;

	  Weights = Vector(3);
	  Xi = Vector(3);
	  Eta = Vector(3);

	  for(int i=0; i<Points; ++i)
	    Weights[i] = 1./3;
	  
	  Xi[0] = .5; Eta[0] = 0.;
	  Xi[1] = .5; Eta[1] = .5;
	  Xi[2] = 0.; Eta[2] = .5;

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
      return Points;
    }
    
    void QuadratureFormula::FormulaData(Vector& weights, Vector& xi, Vector& eta)
    {
      weights = Weights;
      xi = Xi;
      eta = Eta;
    }

  }
}
