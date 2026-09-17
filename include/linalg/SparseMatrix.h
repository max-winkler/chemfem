#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#include <iostream>

#include "linalg/Vector.h"

namespace chemfem{
  namespace linalg{

    enum LIN_SOLVER {CG, GMRES, UMFPACK};

    /**
     * Data structure for sparse matrices stored in compressed row format.
     */
    class SparseMatrix
    {
      friend class IdentityMatrix;
      friend class SparseMatrixInserter;
      friend class DirectSolver;

    public:
      /**
       * Creates an empty matrix of dimension m-n
       */
      SparseMatrix(size_t, size_t);

      /**
       * Generates a hard copy of a matrix
       */
      SparseMatrix(const SparseMatrix&);

      /**
       * Frees the three arrays of the compressed row format. Without it every assembly of a
       * form leaks its whole matrix, which nobody noticed as long as a system was assembled
       * once, and which a Newton iteration turns into gigabytes.
       */
      ~SparseMatrix();

      /**
       * Copies a matrix over an existing one, freeing what it held before
       */
      SparseMatrix& operator=(const SparseMatrix&);

      /**
       * Takes the arrays over from a temporary instead of copying them. This is the case of
       * Matrix = SparseMatrix(n, n), which every assembly of a block system starts with.
       */
      SparseMatrix(SparseMatrix&&) noexcept;
      SparseMatrix& operator=(SparseMatrix&&) noexcept;

      /**
       * Used to print the matrix to the console or write into a file.
       */
      friend std::ostream& operator<<(std::ostream&, const SparseMatrix&);

      /**
       * Multiplication with a vector.
       */
      Vector operator*(const Vector&) const;
      
      /**
       * Solves the linear equation system. The solver is picked from the enumeration
       * LIN_SOLVER. CG needs the matrix to be symmetric and positive definite, GMRES
       * and UMFPACK do not. UMFPACK is direct, so it returns the exact solution up to
       * round off, at the price of storing the factors.
       */
      Vector Solve(const Vector&, LIN_SOLVER = CG);

      /**
       * Sets the entries of a row and a column to zero and the diagonal entry to one, which
       * fixes the corresponding unknown. The diagonal entry has to be part of the sparsity
       * pattern already.
       */
      void EliminateRowAndColumn(size_t);
      
    private:
      size_t *Col = NULL, *Row = NULL;
      double *Entry = NULL;
      size_t m=0, n=0, nnz=0;

      /**
       * Solve equation system with a CG method.
       */
      void Solve_CG(const Vector&, Vector&);

      /**
       * Solve equation system with a restarted GMRES method. Unlike CG this does not
       * rely on symmetry, so it also works for the convection term.
       */
      void Solve_GMRES(const Vector&, Vector&);

      /**
       * Solve equation system with the direct solver UMFPACK.
       */
      void Solve_UMFPACK(const Vector&, Vector&);
    };
    
    /**
     * Returns an Identity Matrix of size n-n
     */
    class IdentityMatrix : public SparseMatrix{
    public:
      /**
       * Initializes an identity matrix with dimension n-n
       */
      IdentityMatrix(size_t);
    };
    
    /**
     * Prints the matrix to the console
     */
    std::ostream& operator<<(std::ostream&, const SparseMatrix&);

  };
};



#endif
