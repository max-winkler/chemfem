#ifndef _DENSE_MATRIX_H_
#define _DENSE_MATRIX_H_

#include <iostream>

namespace chemfem{
  namespace mesh{
    class Cell;
  }
}

namespace chemfem{
  namespace linalg{

    class DenseMatrix
    {
      friend class chemfem::mesh::Cell;
      
    public:
      DenseMatrix(size_t, size_t);
      DenseMatrix(const DenseMatrix&);
      DenseMatrix& operator=(DenseMatrix&);
      ~DenseMatrix();

      double Determinant();
      DenseMatrix Invert();
      
    private:
      size_t m, n;
      double* data;

      void copy(const DenseMatrix&);
    };
    
  };
};

#endif
