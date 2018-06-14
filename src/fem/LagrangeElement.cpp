#include "fem/LagrangeElement.h"

namespace chemfem{
  namespace fem{
    
    LagrangeElement::LagrangeElement(int degree) : Element(FEType::Lagrange, degree)
    {
      switch(degree)
	{
	case 1: nr_dof = 3; break;
	case 2: nr_dof = 6; break;
	case 3: nr_dof = 10; break;
	default:
	  std::cerr << "Lagrange elements of order " << degree << " are not implemented yet\n";  
	}
    }
    
    double LagrangeElement::Value(int i, double x, double y) const
    {
      double lambda[3] = {1.-x-y, x, y};
      switch(degree)
	{
	case 1:
	  if(i<3)
	    return lambda[i];
	  break;
	case 2:
	  if(i<3)
	    return lambda[i]*(2*lambda[i]-1);
	  else if(i<6)
	    return 4*lambda[i-3]*lambda[(i-2)%3];
	  break;
	case 3:
	  if(i<3)
	    return 4.5*lambda[i]*(lambda[i]-1./3)*(lambda[i]-2./3);
	  else if(i<9)
	    {
	      int edge_ind = (i-3) / 2;
	      int vert_ind = (i-3) % 2;  
	      return 13.5 * lambda[i]*lambda[(i+1)%3]*(lambda[(edge_ind+vert_ind)%3]-1./3);
	    }
	  else if(i==9)
	    return 27*lambda[0]*lambda[1]*lambda[2];
	  break;

	}
      std::cerr << "Requested function value of trial function " << i
		<< ", but the element has only " << nr_dof << " degrees of freedom.\n";

      return 0.;
    }

    Vector LagrangeElement::Gradient(int i, double x, double y) const
    {           
      double lambda[3] = {1.-x-y, x, y};
      Vector gradLambda[3] = {Vector(2), Vector(2), Vector(2)};

      Vector grad(2);
      
      gradLambda[0][0] = -1.;
      gradLambda[0][1] = -1.;

      gradLambda[1][0] = 1.;
      gradLambda[1][1] = 0.;

      gradLambda[2][0] = 0.;
      gradLambda[2][1] = 1.;
      
      switch(degree)
	{
	case 1:
	  grad = gradLambda[i];
	  return grad;
	case 2:
	  if(i<3)
	    return (2*lambda[i]-1.)*gradLambda[i] + 2*lambda[i]*gradLambda[i];
	  else
	    return 4*(lambda[i-3]*gradLambda[(i-2)%3] + lambda[(i-2)%3]*gradLambda[i-3]);
	}

      return Vector();
    }
    
  }
}
