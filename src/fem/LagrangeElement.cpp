#include "fem/LagrangeElement.h"

#include <cmath>

#include "linalg/DenseMatrix.h"

using chemfem::linalg::DenseMatrix;

namespace chemfem{
  namespace fem{
    
    LagrangeElement::LagrangeElement(int degree) : Element(FEType::Lagrange, degree)
    {
      switch(degree)
	{
	case 1: nr_dof = 3; break;
	case 2: nr_dof = 6; break;
	case 3: nr_dof = 10; break;
	case 4: nr_dof = 15; break;
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
	      return 13.5 * lambda[edge_ind]*lambda[(edge_ind+1)%3]
		*(lambda[(edge_ind+vert_ind)%3]-1./3);
	    }
	  else if(i==9)
	    return 27*lambda[0]*lambda[1]*lambda[2];
	  break;
	  
	case 4:
	  if(i<3)
	    return 32./3*lambda[i]*(lambda[i]-0.25)*(lambda[i]-0.5)*(lambda[i]-0.75);
	  else if(i<12)
	    {
	      int edge_ind = (i-3) / 3;
	      int vert_ind = (i-3) % 3;
	      if(vert_ind == 0)
		return 128./3*lambda[edge_ind]*(lambda[edge_ind]-0.25)
		  * (lambda[edge_ind]-0.5)*lambda[(edge_ind+1)%3];
	      else if(vert_ind == 1)
		return 64*lambda[edge_ind]*(lambda[edge_ind]-0.25)
		  * lambda[(edge_ind+1)%3]*(lambda[(edge_ind+1)%3]-0.25);
	      if(vert_ind == 2)
		return 128./3*lambda[(edge_ind+1)%3]*(lambda[(edge_ind+1)%3]-0.25)
		  * (lambda[(edge_ind+1)%3]-0.5)*lambda[edge_ind];
	    }
	  else if(i<15)
	    {
	      int int_ind = i-12;
	      return 128.*lambda[0]*lambda[1]*lambda[2]*(lambda[int_ind]-0.25);
	    }
	  break;

	}
      std::cerr << "Requested function value of trial function " << i
		<< ", but the element has only " << nr_dof << " degrees of freedom.\n";

      return 0.;
    }

    Vector LagrangeElement::Gradient(int i, double x, double y) const
    {           
      double lambda[3] = {1.-x-y, x, y};

      if(x > 1 || y > 1 || x < 0 || y < 0 || x+y > 1)
	std::cerr << "Invalid quadrature point.\n";
      
      Vector grad_L(3);
            
      switch(degree)
	{
	case 1:
	  if(i<3)
	    grad_L[i] = 1.;
	  break;
	  
	case 2:
	  if(i<3)
	    grad_L[i] = 4.*lambda[i] -1.;
	  else if(i<6)
	    {
	      grad_L[i-3] = 4.*lambda[(i-2)%3];
	      grad_L[(i-2)%3] = 4.*lambda[i-3];
	    }
	  break;
	  
	case 3:
	  if(i<3)
	    grad_L[i] = 13.5*lambda[i]*lambda[i] - 9*lambda[i] + 1.;
	  else if(i<9)
	    {
	      int edge_ind = (i-3) / 2;
	      int vert_ind = (i-3) % 2;

	      grad_L[edge_ind] = 13.5*lambda[(edge_ind+1)%3]
		*(lambda[(edge_ind + vert_ind)%3]-1./3 + (vert_ind == 0 ? lambda[edge_ind] : 0.));
	      grad_L[(edge_ind+1)%3] = 13.5*lambda[edge_ind]
		*(lambda[(edge_ind + vert_ind)%3]-1./3 + (vert_ind == 1 ? lambda[(edge_ind+1)%3] : 0.));
	    }
	  else if(i==9)
	    {
	      grad_L[0] = 27.*lambda[1]*lambda[2];
	      grad_L[1] = 27.*lambda[0]*lambda[2];
	      grad_L[2] = 27.*lambda[0]*lambda[1];
	    }
	  break;

	case 4:
	  if(i<3)
	    grad_L[i] = 128./3*pow(lambda[i], 3) - 48.*pow(lambda[i], 2) + 44./3*lambda[i] - 1.;
	  else if(i<12)
	    {
	      int edge_ind = (i-3) / 3;
	      int vert_ind = (i-3) % 3;
	      if(vert_ind == 0)
		{
		  grad_L[edge_ind] = 128./3*lambda[(edge_ind+1)%3]
		    * (3*lambda[edge_ind]*lambda[edge_ind] - 1.5*lambda[edge_ind] + 1./8);
		  grad_L[(edge_ind+1)%3] = 128./3*lambda[edge_ind]
		    * (lambda[edge_ind]-0.25)*(lambda[edge_ind]-0.5);		  
		}
	      else if(vert_ind == 1)
		{
		  grad_L[edge_ind] = 64.*lambda[(edge_ind+1)%3]*(lambda[(edge_ind+1)%3]-0.25)
		    * (2*lambda[edge_ind]-0.25);
		  grad_L[(edge_ind+1)%3] = 64.*lambda[edge_ind]*(lambda[edge_ind]-0.25)
		    * (2*lambda[(edge_ind+1)%3]-0.25);
		}
	      else if(vert_ind == 2)
		{
		  grad_L[edge_ind] = 128./3*lambda[(edge_ind+1)%3]
		    * (lambda[(edge_ind+1)%3]-0.25)*(lambda[(edge_ind+1)%3]-0.5);
		  grad_L[(edge_ind+1)%3] = 128./3*lambda[edge_ind]
		    * (3*lambda[(edge_ind+1)%3]*lambda[(edge_ind+1)%3]
		       - 1.5*lambda[(edge_ind+1)%3] + 1./8);
		}
	    }
	  else if(i<15)
	    {
	      int int_ind = i-12;
	      grad_L[int_ind] = 128.*lambda[(int_ind+1)%3]*lambda[(int_ind+2)%3]
		* (2*lambda[int_ind] - 0.25);
	      grad_L[(int_ind+1)%3] = 128.*lambda[(int_ind+2)%3]*lambda[int_ind]
		* (lambda[int_ind]-0.25);
	      grad_L[(int_ind+2)%3] = 128.*lambda[(int_ind+1)%3]*lambda[int_ind]
		* (lambda[int_ind]-0.25);
	    }
	  break;
	}

      // Chain rule with the derivative of the barycentric coordinates,
      //   d(lambda_0,lambda_1,lambda_2)/d(x,y) = ( -1  1  0 )
      //                                          ( -1  0  1 ).
      // That matrix is constant, so the product is written out here instead of
      // building a DenseMatrix on every call.
      Vector grad(2);
      grad[0] = grad_L[1] - grad_L[0];
      grad[1] = grad_L[2] - grad_L[0];

      return grad;
    }
    
  }
}
