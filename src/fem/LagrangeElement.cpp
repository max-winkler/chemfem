#include "fem/LagrangeElement.h"

namespace chemfem{
  namespace fem{
    
    LagrangeElement::LagrangeElement(int degree) : Element(FEType::Lagrange, degree)
    {
      if(degree > 1)
	{
	  std::cerr << "Lagrange elements of order " << degree << " are not implemented yet\n";  
	}
    }
    
    double LagrangeElement::Value(int i, double x, double y) const
    {
      switch(degree)
	{
	case 1:
	  switch(i)
	    {
	    case 0: return 1.-x-y;
	    case 1: return x;
	    case 2: return y;
	    }
	  break;
	}
    }

    Vector LagrangeElement::Gradient(int i, double x, double y) const
    {
      Vector grad(2);
      switch(degree)
	{
	case 1:
	  switch(i)
	    {
	    case 0: grad[0] = -1.; grad[1] = -1.; break;
	    case 1: grad[0] = 1.; grad[1] = 0.; break;
	    case 2: grad[0] = 0.; grad[1] = 1.; break;
	    }
	  return grad;	  
	}      
    }    
    
  }
}
