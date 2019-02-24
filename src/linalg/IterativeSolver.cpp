#include <cmath>

#include "linalg/SparseMatrix.h"

namespace chemfem{
  namespace linalg{

    Vector SparseMatrix::Solve(const Vector& b)
    {
      if(m != n)
	{
	  std::cerr << "Matrix is not quadratic. Cannot solve equation system.\n";
	  return Vector();
	}

      Vector x(m);
      Solve_CG(b, x);

      return x;
    }

    void SparseMatrix::Solve_CG(const Vector& b, Vector& x)
    {
      Vector res(m);

      res = b - (*this)*x;

      Vector d(res);

      const int max_iter = 1000;
      const double a_tol = 1.e-8;

      int iter = 0;

      double eps, eps_old;
      
      while(iter < max_iter)
	{
	  Vector z((*this)*d);
	  eps = dot(res, res);

	  if(iter % 100 == 0)
	    {
	      std::cout << "CG it " << iter << ", resid = " << sqrt(eps) << "                ";
	      std::cout.flush();
	      std::cout << "\r";
	    }
	  
	  double alpha = eps / dot(d, z);
	  d.axpy(alpha, x, x);
	  z.axpy(-alpha, res, res);

	  eps_old = eps;
	  eps = dot(res, res);
	  double beta = eps / eps_old;

	  d.axpy(beta, res, d);

	  if(sqrt(eps) < a_tol)
	    break;

	  ++iter;
	}
      std::cout << "CG finished after " << iter << " iterations with resid = " << sqrt(eps) << std::endl;
    }
        
  }
}
