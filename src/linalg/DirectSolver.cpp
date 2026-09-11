#include <climits>
#include <iostream>
#include <vector>

#include <umfpack.h>

#include "linalg/DirectSolver.h"

namespace chemfem{
  namespace linalg{

    // UMFPACK stores a matrix by columns and indexes it with int, SparseMatrix stores it by
    // rows and indexes it with size_t. The compressed row storage of A is however exactly
    // the compressed column storage of A^T, so the arrays are handed over unchanged and
    // UMFPACK factorizes A^T. The transposed solve with UMFPACK_At then yields the solution
    // of A x = b again. The flag is therefore not optional: with UMFPACK_A the result would
    // silently be wrong for every non symmetric matrix.

    DirectSolver::DirectSolver(const SparseMatrix& A) : n(A.m), Numeric(NULL)
    {
      if(A.m != A.n)
	{
	  std::cerr << "Error: UMFPACK needs a square matrix.\n";
	  return;
	}

      if(A.m+1 > (size_t)INT_MAX || A.nnz > (size_t)INT_MAX)
	{
	  std::cerr << "Error: The matrix is too large for the int indices of UMFPACK.\n";
	  return;
	}

      Ap.assign(A.Row, A.Row + A.m + 1);
      Ai.assign(A.Col, A.Col + A.nnz);
      Ax.assign(A.Entry, A.Entry + A.nnz);

      void *Symbolic = NULL;
      double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];

      umfpack_di_defaults(Control);

      int status = umfpack_di_symbolic((int)A.n, (int)A.m, Ap.data(), Ai.data(), Ax.data(),
				       &Symbolic, Control, Info);

      if(status != UMFPACK_OK)
	{
	  std::cerr << "Error: The symbolic factorization of UMFPACK failed with status "
		    << status << ".\n";
	  umfpack_di_free_symbolic(&Symbolic);
	  return;
	}

      status = umfpack_di_numeric(Ap.data(), Ai.data(), Ax.data(), Symbolic, &Numeric,
				  Control, Info);

      umfpack_di_free_symbolic(&Symbolic);

      if(status != UMFPACK_OK)
	{
	  std::cerr << "Error: The numeric factorization of UMFPACK failed with status "
		    << status;
	  if(status == UMFPACK_WARNING_singular_matrix)
	    std::cerr << ", the matrix is singular";
	  std::cerr << ".\n";

	  umfpack_di_free_numeric(&Numeric);
	  Numeric = NULL;
	}
    }

    DirectSolver::~DirectSolver()
    {
      if(Numeric)
	umfpack_di_free_numeric(&Numeric);
    }

    Vector DirectSolver::Solve(const Vector& b) const
    {
      Vector x(n);

      if(!Numeric)
	{
	  std::cerr << "Error: There is no valid factorization to solve with.\n";
	  return x;
	}

      if(b.size() != n)
	{
	  std::cerr << "Error: The right hand side has length " << b.size()
		    << ", but the matrix has " << n << " rows.\n";
	  return x;
	}

      double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
      umfpack_di_defaults(Control);

      int status = umfpack_di_solve(UMFPACK_At, Ap.data(), Ai.data(), Ax.data(), &x[0], &b[0],
				    Numeric, Control, Info);

      if(status != UMFPACK_OK)
	std::cerr << "Error: The solve step of UMFPACK failed with status " << status << ".\n";

      return x;
    }

    void SparseMatrix::Solve_UMFPACK(const Vector& b, Vector& x)
    {
      DirectSolver LU(*this);
      x = LU.Solve(b);
    }

  }
}
