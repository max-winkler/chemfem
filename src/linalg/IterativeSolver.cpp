#include <algorithm>
#include <cmath>
#include <vector>

#include "linalg/SparseMatrix.h"

namespace chemfem{
  namespace linalg{

    Vector SparseMatrix::Solve(const Vector& b, LIN_SOLVER solver)
    {
      if(m != n)
	{
	  std::cerr << "Matrix is not quadratic. Cannot solve equation system.\n";
	  return Vector();
	}

      if(b.size() != m)
	{
	  std::cerr << "Right hand side has length " << b.size() << ", but the matrix has "
		    << m << " rows.\n";
	  return Vector();
	}

      Vector x(m);

      switch(solver)
	{
	case CG:
	  Solve_CG(b, x);
	  break;
	case GMRES:
	  Solve_GMRES(b, x);
	  break;
	case UMFPACK:
	  Solve_UMFPACK(b, x);
	  break;
	}

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

    void SparseMatrix::Solve_GMRES(const Vector& b, Vector& x)
    {
      // Length of a Krylov cycle. Everything is restarted after this many steps, which
      // caps the storage at restart+1 vectors and the orthogonalization at restart
      // inner products per step.
      const size_t restart = 30;

      const int max_iter = 1000;
      const double a_tol = 1.e-8;

      std::vector<Vector> V(restart+1, Vector(m));

      // Hessenberg matrix of the Arnoldi process, (restart+1) x restart, stored by rows
      std::vector<double> H((restart+1)*restart, 0.);

      // Givens rotations that keep it upper triangular, and the rotated right hand side
      std::vector<double> cs(restart, 0.), sn(restart, 0.), g(restart+1, 0.);

      int iter = 0;
      double resid = 0.;

      while(iter < max_iter)
	{
	  Vector r(b - (*this)*x);

	  const double beta = sqrt(dot(r, r));
	  resid = beta;

	  if(beta < a_tol)
	    break;

	  V[0] = (1./beta) * r;

	  std::fill(g.begin(), g.end(), 0.);
	  g[0] = beta;

	  size_t k = 0;

	  for(size_t j=0; j<restart && iter<max_iter; ++j, ++iter)
	    {
	      Vector w((*this)*V[j]);

	      // Arnoldi step, orthogonalized against the basis built so far. Modified
	      // Gram-Schmidt subtracts one direction at a time, which loses far less
	      // orthogonality in finite precision than the classical variant.
	      for(size_t i=0; i<=j; ++i)
		{
		  H[i*restart + j] = dot(w, V[i]);
		  V[i].axpy(-H[i*restart + j], w, w);
		}

	      const double h_next = sqrt(dot(w, w));
	      H[(j+1)*restart + j] = h_next;

	      if(h_next > 0.)
		V[j+1] = (1./h_next) * w;

	      // Apply the rotations of the previous steps to the new column
	      for(size_t i=0; i<j; ++i)
		{
		  const double t = H[i*restart + j];

		  H[i*restart + j]     =  cs[i]*t + sn[i]*H[(i+1)*restart + j];
		  H[(i+1)*restart + j] = -sn[i]*t + cs[i]*H[(i+1)*restart + j];
		}

	      // New rotation that eliminates the subdiagonal entry
	      const double d = sqrt(H[j*restart + j]*H[j*restart + j] + h_next*h_next);

	      if(d == 0.)
		break;

	      cs[j] = H[j*restart + j] / d;
	      sn[j] = h_next / d;

	      H[j*restart + j] = d;
	      H[(j+1)*restart + j] = 0.;

	      // The same rotation on the right hand side. Its last entry is the residual
	      // norm, so it comes for free without forming the iterate.
	      g[j+1] = -sn[j]*g[j];
	      g[j]   =  cs[j]*g[j];

	      resid = fabs(g[j+1]);
	      k = j+1;

	      if(iter % 100 == 0)
		{
		  std::cout << "GMRES it " << iter << ", resid = " << resid << "                ";
		  std::cout.flush();
		  std::cout << "\r";
		}

	      if(resid < a_tol)
		break;
	    }

	  if(k == 0)
	    break;

	  // Back substitution for the triangular system, the result overwrites g
	  for(size_t i=k; i-- > 0; )
	    {
	      double sum = g[i];

	      for(size_t l=i+1; l<k; ++l)
		sum -= H[i*restart + l]*g[l];

	      g[i] = sum / H[i*restart + i];
	    }

	  for(size_t i=0; i<k; ++i)
	    V[i].axpy(g[i], x, x);

	  if(resid < a_tol)
	    break;
	}

      std::cout << "GMRES finished after " << iter << " iterations with resid = "
		<< resid << std::endl;
    }
        
  }
}
