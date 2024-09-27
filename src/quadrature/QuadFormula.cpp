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
	    Weights[i] = 1./6;
	  
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

	case GAUSS_7:
	  Points = 16;

	  Weights = Vector(16);
	  Xi = Vector(16);
	  Eta = Vector(16);

	  Xi[0] = 0.0571041961; Eta[0] = 0.06546699455602246;
	  Xi[1] = 0.2768430136; Eta[1] = 0.05021012321401679;
	  Xi[2] = 0.5835904324; Eta[2] = 0.02891208422223085;
	  Xi[3] = 0.8602401357; Eta[3] = 0.009703785123906346;
	  Xi[4] = 0.0571041961; Eta[4] = 0.3111645522491480;   
	  Xi[5] = 0.2768430136; Eta[5] = 0.2386486597440242;    
	  Xi[6] = 0.5835904324; Eta[6] = 0.1374191041243166;    
	  Xi[7] = 0.8602401357; Eta[7] = 0.04612207989200404;
	  Xi[8] = 0.0571041961; Eta[8] = 0.6317312516508520;    
	  Xi[9] = 0.2768430136; Eta[9] = 0.4845083266559759;    
	  Xi[10] = 0.5835904324; Eta[10] = 0.2789904634756834;    
	  Xi[11] = 0.8602401357; Eta[11] = 0.09363778440799593;
	  Xi[12] = 0.0571041961; Eta[12] = 0.8774288093439775;    
	  Xi[13] = 0.2768430136; Eta[13] = 0.6729468631859832;    
	  Xi[14] = 0.5835904324; Eta[14] = 0.3874974833777692;    
	  Xi[15] = 0.8602401357; Eta[15] = 0.1300560791760936;    

	  Weights[0] = 0.04713673637581137;
	  Weights[1] = 0.07077613579259895;
	  Weights[2] = 0.04516809856187617;
	  Weights[3] = 0.01084645180365496;
	  Weights[4] = 0.08837017702418863;
	  Weights[5] = 0.1326884322074010;
	  Weights[6] = 0.08467944903812383;
	  Weights[7] = 0.02033451909634504;
	  Weights[8] = 0.08837017702418863;
	  Weights[9] = 0.1326884322074010;
	  Weights[10] = 0.08467944903812383;
	  Weights[11] = 0.02033451909634504;
	  Weights[12] = 0.04713673637581137;
	  Weights[13] = 0.07077613579259895;
	  Weights[14] = 0.04516809856187617;
	  Weights[15] = 0.01084645180365496;

	  break;

    case LINE_GAUSS_5:

      Weights[0] = 0.5 * 0.236926885;
      Weights[1] = 0.5 * 0.478628670;
      Weights[2] = 0.5 * 0.568888889;
      Weights[3] = 0.5 * 0.478628670;
      Weights[4] = 0.5 * 0.236926885;

      Xi[0] = 0.5 * (1 - 0.906179846);
      Xi[1] = 0.5 * (1 - 0.538469310);
      Xi[2] = 0.5 * (1 + 0.000000000);
      Xi[3] = 0.5 * (1 + 0.538469310);
      Xi[3] = 0.5 * (1 + 0.906179846);
      
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
