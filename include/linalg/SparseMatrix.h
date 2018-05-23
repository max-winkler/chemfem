#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#include <iostream>

namespace chemfem{
  namespace linalg{

    /**
     * Data structure for sparse matrices stored in compressed row format.
     */
    class SparseMatrix
    {
      friend class IdentityMatrix;
      friend class SparseMatrixInserter;
      
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
       *
       */
      friend std::ostream& operator<<(std::ostream&, const SparseMatrix&);

    private:
      size_t *Col = NULL, *Row = NULL;
      double *Entry = NULL;
      size_t m=0, n=0, nnz=0;

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
