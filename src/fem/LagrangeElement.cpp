#include "fem/LagrangeElement.h"

#include <cmath>
#include <iostream>

#include "linalg/Vector2D.h"

using chemfem::linalg::Coordinate;
using chemfem::linalg::Vector2D;

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

      dofs_per_vertex = 1;
      dofs_per_edge = degree-1;
      dofs_interior = (degree-1)*(degree-2)/2;
    }

    Coordinate LagrangeElement::NodalPoint(int i) const
    {
      const Coordinate Vertex[3] = {{0., 0.}, {1., 0.}, {0., 1.}};

      if(i < 3)
	return Vertex[i];

      if(i < 3 + 3*dofs_per_edge)
	{
	  const int edge = (i-3) / dofs_per_edge;
	  const int j = (i-3) % dofs_per_edge;

	  double xi, eta;
	  EdgeToRefCoords(edge, double(j+1)/degree, xi, eta);

	  return Coordinate{xi, eta};
	}

      const int j = i - 3 - 3*dofs_per_edge;

      switch(degree)
	{
	case 3:
	  return Coordinate{1./3, 1./3};

	case 4:
	  {
	    // Barycentric coordinate 1/2 for the vertex j, 1/4 for the other two
	    double lambda[3] = {0.25, 0.25, 0.25};
	    lambda[j] = 0.5;
	    return Coordinate{lambda[1], lambda[2]};
	  }
	}

      std::cerr << "Requested nodal point " << i << ", but the element has only "
		<< nr_dof << " degrees of freedom.\n";

      return Coordinate{0., 0.};
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

    Vector2D LagrangeElement::Gradient(int i, double x, double y) const
    {           
      double lambda[3] = {1.-x-y, x, y};

      if(x > 1 || y > 1 || x < 0 || y < 0 || x+y > 1)
	std::cerr << "Invalid quadrature point.\n";
      
      // Derivatives with respect to the three barycentric coordinates
      double grad_L[3] = {0., 0., 0.};

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

      // Chain rule with  d(lambda_0,lambda_1,lambda_2)/d(x,y) = ( -1  1  0 )
      //                                                           ( -1  0  1 )
      return Vector2D(grad_L[1] - grad_L[0], grad_L[2] - grad_L[0]);
    }

    Matrix2D LagrangeElement::Hessian(int i, double x, double y) const
    {
      double lambda[3] = {1.-x-y, x, y};

      if(x > 1 || y > 1 || x < 0 || y < 0 || x+y > 1)
	std::cerr << "Invalid quadrature point.\n";

      // Second derivatives with respect to the three barycentric coordinates
      double hess_L[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

      switch(degree)
	{
	case 1:
	  // Affine in the barycentric coordinates, so the second derivatives vanish
	  break;

	case 2:
	  if(i<3)
	    hess_L[i][i] = 4.;
	  else if(i<6)
	    {
	      int a = i-3;
	      int b = (i-2)%3;
	      hess_L[a][b] = hess_L[b][a] = 4.;
	    }
	  break;

	case 3:
	  if(i<3)
	    hess_L[i][i] = 27.*lambda[i] - 9.;
	  else if(i<9)
	    {
	      int edge_ind = (i-3) / 2;
	      int vert_ind = (i-3) % 2;

	      // The quadratic factor sits on the vertex p of the edge, the remaining
	      // vertex o of that edge enters linearly
	      int p = (edge_ind + vert_ind) % 3;
	      int o = (edge_ind + 1 - vert_ind) % 3;

	      hess_L[p][p] = 27.*lambda[o];
	      hess_L[p][o] = hess_L[o][p] = 13.5*(2.*lambda[p] - 1./3);
	    }
	  else if(i==9)
	    {
	      hess_L[0][1] = hess_L[1][0] = 27.*lambda[2];
	      hess_L[0][2] = hess_L[2][0] = 27.*lambda[1];
	      hess_L[1][2] = hess_L[2][1] = 27.*lambda[0];
	    }
	  break;

	case 4:
	  if(i<3)
	    hess_L[i][i] = 128.*lambda[i]*lambda[i] - 96.*lambda[i] + 44./3;
	  else if(i<12)
	    {
	      int edge_ind = (i-3) / 3;
	      int vert_ind = (i-3) % 3;
	      int next_ind = (edge_ind+1) % 3;

	      if(vert_ind == 0)
		{
		  hess_L[edge_ind][edge_ind] = 128./3*lambda[next_ind]
		    * (6.*lambda[edge_ind] - 1.5);
		  hess_L[edge_ind][next_ind] = hess_L[next_ind][edge_ind] = 128./3
		    * (3*lambda[edge_ind]*lambda[edge_ind] - 1.5*lambda[edge_ind] + 1./8);
		}
	      else if(vert_ind == 1)
		{
		  hess_L[edge_ind][edge_ind] = 128.*lambda[next_ind]*(lambda[next_ind]-0.25);
		  hess_L[next_ind][next_ind] = 128.*lambda[edge_ind]*(lambda[edge_ind]-0.25);
		  hess_L[edge_ind][next_ind] = hess_L[next_ind][edge_ind]
		    = 64.*(2*lambda[edge_ind]-0.25)*(2*lambda[next_ind]-0.25);
		}
	      else if(vert_ind == 2)
		{
		  hess_L[next_ind][next_ind] = 128./3*lambda[edge_ind]
		    * (6.*lambda[next_ind] - 1.5);
		  hess_L[edge_ind][next_ind] = hess_L[next_ind][edge_ind] = 128./3
		    * (3*lambda[next_ind]*lambda[next_ind] - 1.5*lambda[next_ind] + 1./8);
		}
	    }
	  else if(i<15)
	    {
	      int int_ind = i-12;
	      int next_ind = (int_ind+1) % 3;
	      int last_ind = (int_ind+2) % 3;

	      hess_L[int_ind][int_ind] = 256.*lambda[next_ind]*lambda[last_ind];
	      hess_L[int_ind][next_ind] = hess_L[next_ind][int_ind]
		= 128.*lambda[last_ind]*(2*lambda[int_ind] - 0.25);
	      hess_L[int_ind][last_ind] = hess_L[last_ind][int_ind]
		= 128.*lambda[next_ind]*(2*lambda[int_ind] - 0.25);
	      hess_L[next_ind][last_ind] = hess_L[last_ind][next_ind]
		= 128.*lambda[int_ind]*(lambda[int_ind] - 0.25);
	    }
	  break;
	}

      // Chain rule with  d(lambda_0,lambda_1,lambda_2)/d(x,y) = ( -1  1  0 )
      //                                                         ( -1  0  1 )
      // The barycentric coordinates are affine in (x,y), so their own second
      // derivatives drop out and only the congruence with that matrix remains.
      const double grad_lambda[3][2] = {{-1., -1.}, {1., 0.}, {0., 1.}};

      double H[2][2] = {{0., 0.}, {0., 0.}};

      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  for(int k=0; k<2; ++k)
	    for(int l=0; l<2; ++l)
	      H[k][l] += hess_L[a][b] * grad_lambda[a][k] * grad_lambda[b][l];

      return Matrix2D(H[0][0], H[0][1], H[1][0], H[1][1]);
    }
    
  }
}
