#include "fem/LagrangeElement.h"

namespace chemfem{
  namespace fem{
    
    LagrangeElement::LagrangeElement(int degree) : Element(FEType::Lagrange, degree)
    {}
    
    double LagrangeElement::Value(int i, double x, double y)
    {
      switch(i)
	{
	case 0:
	  return 1.-x-y;
	case 1:
	  return x;
	case 2:
	  return y;
	}
    }
  }
}
