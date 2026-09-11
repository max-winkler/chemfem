#include <algorithm>
#include <climits>
#include <vector>

#include <umfpack.h>

#include "linalg/SparseMatrix.h"

namespace chemfem{
  namespace linalg{

    void SparseMatrix::Solve_UMFPACK(const Vector& b, Vector& x)
    {
      // UMFPACK stores a matrix by columns and indexes it with int, this class stores
      // it by rows and indexes it with size_t. The compressed row storage of A is
      // however exactly the compressed column storage of A^T, so the arrays are handed
      // over unchanged and UMFPACK factorizes A^T. Asking it for the transposed solve
      // with UMFPACK_At then yields the solution of A x = b again. That saves building
      // an explicit transpose, but it does mean the flag is not optional: with
      // UMFPACK_A the result would silently be wrong for every non symmetric matrix.
      if(m+1 > (size_t)INT_MAX || nnz > (size_t)INT_MAX)
	{
	  std::cerr << "Error: The matrix is too large for the int indices of UMFPACK.\n";
	  return;
	}

      std::vector<int> Ap(m+1), Ai(nnz);

      std::copy(Row, Row+m+1, Ap.begin());
      std::copy(Col, Col+nnz, Ai.begin());

      void *Symbolic = NULL, *Numeric = NULL;
      double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];

      umfpack_di_defaults(Control);

      int status = umfpack_di_symbolic((int)n, (int)m, &Ap[0], &Ai[0], Entry,
				       &Symbolic, Control, Info);

      if(status != UMFPACK_OK)
	{
	  std::cerr << "Error: The symbolic factorization of UMFPACK failed with status "
		    << status << ".\n";
	  umfpack_di_free_symbolic(&Symbolic);
	  return;
	}

      status = umfpack_di_numeric(&Ap[0], &Ai[0], Entry, Symbolic, &Numeric, Control, Info);

      umfpack_di_free_symbolic(&Symbolic);

      if(status != UMFPACK_OK)
	{
	  std::cerr << "Error: The numeric factorization of UMFPACK failed with status "
		    << status;
	  if(status == UMFPACK_WARNING_singular_matrix)
	    std::cerr << ", the matrix is singular";
	  std::cerr << ".\n";

	  umfpack_di_free_numeric(&Numeric);
	  return;
	}

      status = umfpack_di_solve(UMFPACK_At, &Ap[0], &Ai[0], Entry, &x[0], &b[0],
				Numeric, Control, Info);

      umfpack_di_free_numeric(&Numeric);

      if(status != UMFPACK_OK)
	std::cerr << "Error: The solve step of UMFPACK failed with status " << status << ".\n";
    }

  }
}
