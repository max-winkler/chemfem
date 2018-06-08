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
    }

    void SparseMatrix::Solve_CG(const Vector& b, Vector& x)
    {
      Vector resid(m);

      // resid = b - (*this)*x;
    }
    
  }
}
