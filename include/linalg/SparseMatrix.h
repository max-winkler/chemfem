#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

namespace chemfem{
  namespace linalg{

    /**
     * Data structure for sparse matrices stored in compressed row format.
     */
    class SparseMatrix
    {
    private:
      size_t *Col, *Row;
      double *Entry;
      size_t m, n;
      
    public:
      /**
       * Creates an empty matrix of dimension m-n
       */
      void SparseMatrix(size_t, size_t);

      /**
       * Generates a hard copy of a matrix
       */
      void SparseMatrix(const SparseMatrix&);
      
    };
  }
}

#endif
