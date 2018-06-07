#ifndef _DENSE_MATRIX_H_
#define _DENSE_MATRIX_H_

#include <iostream>

#include "linalg/Vector.h"

namespace chemfem{
  namespace mesh{
    class Cell;
  }
}

namespace chemfem{
  namespace linalg{

    /**
     * Each instance of this class represents a dense matrix. The data are stored row-wise.     
     */
    class DenseMatrix
    {
      friend class chemfem::mesh::Cell;

      friend std::ostream& operator<<(std::ostream&, const DenseMatrix&);
      
      class DenseMatrixRow
      {
      public:
	DenseMatrixRow(double*);
	double& operator[](const size_t);
	const double& operator[](const size_t) const;
      private:
	double* RowPtr;
      };
      
    public:
      /**
       * Constructor which initializes the matrix by its dimension.
       */
      DenseMatrix(size_t, size_t);

      /**
       * Copy constructure. The data of the matrix is copied.
       */
      DenseMatrix(const DenseMatrix&);

      /**
       * Copy assignment. The data of the matrix is copied.
       */
      DenseMatrix& operator=(DenseMatrix&);

      /**
       * Frees all allocated memory.
       */
      ~DenseMatrix();

      /**
       * Writes the entries row-wise into the matrix
       */
      void Init(const double*);
      
      /**
       * Returns the determinant of the matrix. Only implemented for 2x2 matrices.
       */
      double Determinant() const;

      /**
       * Returns the inverse of the matrix. Only implemented for 2x2 matrices.
       */
      DenseMatrix Invert() const;

      /**
       * Returns the transposed of the matrix.
       */
      DenseMatrix Transpose() const;
      
      /**
       * Computes a matrix vector product and returns its result.
       */
      Vector operator*(const Vector&) const;

      /**
       * Returns a row of the matrix. DenseMatrixRow is just used to allow the entry of a 
       * single element of the matrix by means of [i][j].
       */
      DenseMatrixRow operator[](const size_t);

      /**
       * Returns a row of the matrix. DenseMatrixRow is just used to allow the entry of a 
       * single element of the matrix by means of [i][j].
       */
      const DenseMatrixRow operator[](const size_t) const;
      
    private:
      size_t m, n;
      double* data;

      void copy(const DenseMatrix&);
    };
  };
};

#endif
