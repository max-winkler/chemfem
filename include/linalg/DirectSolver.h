#ifndef _DIRECT_SOLVER_H_
#define _DIRECT_SOLVER_H_

#include <vector>

#include "linalg/SparseMatrix.h"
#include "linalg/Vector.h"

namespace chemfem{
  namespace linalg{

    /**
     * LU factorization of a sparse matrix with UMFPACK. The matrix is factorized once in
     * the constructor, afterwards each call of Solve() only costs a forward and a backward
     * substitution. This pays off whenever the same matrix is solved with many right hand
     * sides, e.g. in a time stepping scheme with a constant step size.
     *
     * The solver keeps its own copy of the matrix, later changes of the matrix do not
     * affect it.
     */
    class DirectSolver
    {
    public:
      /**
       * Factorizes the matrix. Iterative refinement gains a few digits in each solve and
       * costs about two more triangular solves, so it is worth switching off when many
       * right hand sides are solved with the same matrix.
       */
      DirectSolver(const SparseMatrix&, bool IterativeRefinement = true);

      ~DirectSolver();

      DirectSolver(const DirectSolver&) = delete;
      DirectSolver& operator=(const DirectSolver&) = delete;

      Vector Solve(const Vector&) const;

    private:
      size_t n;
      bool Refinement;
      std::vector<int> Ap, Ai;
      std::vector<double> Ax;
      void* Numeric;
    };

  };
};

#endif
