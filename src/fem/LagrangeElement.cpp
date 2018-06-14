#include "fem/LagrangeElement.h"

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

	}
      std::cerr << "Requested function value of trial function " << i
		<< ", but the element has only " << nr_dof << " degrees of freedom.\n";

      return 0.;
    }

    Vector LagrangeElement::Gradient(int i, double x, double y) const
    {           
      double lambda[3] = {1.-x-y, x, y};

      DenseMatrix dLdX(2,3);
      dLdX[0][0] = -1.; dLdX[0][1] = 1.; dLdX[0][2] = 0.;
      dLdX[1][0] = -1.; dLdX[1][1] = 0.; dLdX[1][2] = 1.;

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
	      grad_L[i] = 4*lambda[(i+1)%3];
	      grad_L[(i+1)%3] = 4*lambda[i];
	    }
	  break;
	case 3:
	  if(i<3)
	    grad_L[i] = 4.5*(lambda[i]*(lambda[i]-1./3) + (2*lambda[i]-1./3)*(lambda[i]-2./3));
	  else if(i<9)
	    {
	      int edge_ind = (i-3) / 2;
	      int vert_ind = (i-3) % 2;

	      grad_L[edge_ind] = 4.5*lambda[(edge_ind+1)%3]
		*(lambda[edge_ind + vert_ind]-1./3 + (vert_ind == 0 ? lambda[edge_ind] : 0.));
	      grad_L[(edge_ind+1)%3] = 4.5*lambda[edge_ind]
		*(lambda[edge_ind + vert_ind]-1./3 + (vert_ind == 1 ? lambda[(edge_ind+1)%3] : 0.));
	    }
	  else if(i==10)
	    {
	      grad_L[0] = 27.*lambda[1]*lambda[2];
	      grad_L[1] = 27.*lambda[0]*lambda[2];
	      grad_L[2] = 27.*lambda[0]*lambda[1];
	    }
	    break;
	}

      return dLdX*grad_L;
    }
    
  }
}
